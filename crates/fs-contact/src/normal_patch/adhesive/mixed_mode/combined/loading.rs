//! First-boundary loading of the original shared I/II/III contact.
//!
//! The apparatus holds normal approach and prescribes
//! q=[delta,lambda*u1,lambda*u2,lambda*phi]. The dimensionless signed lambda
//! scales explicit metre/radian directions; these are not independent patches.
//! With A=E*(a²/R-delta)²-2*pi*W*a and B=4G*|u|²+8Gbar*a²*phi²,
//! equilibrium is A+lambda²*B=0. The original stable opening radius decreases
//! with |lambda| until its first radial fold or normal-opening edge.
//!
//! The scalar owner below brackets that first boundary without scanning
//! branches or inverting the vanishing tangent. Boundary responses are
//! numerical root approximations, with radius brackets and actual original
//! stationarity/curvature residuals exposed. They are not interval certificates.
//! Endpoints and, on reversal, lambda=0 enclose every visited radius. Original
//! normal screens are monotone or convex on that interval; both joint nominal
//! shear/twist screens grow with |lambda|. Thus the complete accepted segment
//! is admitted, including unloading and sign reversal.
//!
//! Work is a factored difference of all four original potential stores. Only
//! an authenticated radial fold at negative approach can release into a
//! detached state. That nonnegative release has an explicit unresolved
//! mechanical destination: no impulse, heat fraction, slip law or rebonding
//! is inferred. Opening-edge events and folds at nonnegative approach stop.

mod radius;

/// Original-source thermal changes with retained loading and release history.
pub mod temperature;

/// Explicit affine apparatus clocks and the first accepted source boundary.
pub mod timed;

use super::{JkrCombinedBoundary, JkrCombinedError, JkrCombinedResponse, SourceJkrMixedMode};
use core::f64::consts::PI;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedLoadingControls {
    /// Scalar root iterations during preparation, with no branch scan.
    pub maximum_iterations: usize,
    /// Actual polynomial/constitutive probes per operation, including failures.
    pub maximum_evaluations: usize,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum JkrCombinedLoadingError {
    Source(JkrCombinedError),
    Invalid(&'static str),
    StateMismatch,
    Budget,
    Cancelled,
    UnresolvedBoundary,
    BeyondFirstBoundary,
    NotAtSeparatingFold,
    ClockRoundoff { residual_s: f64, tolerance_s: f64 },
    EnergyBalance { residual_j: f64, tolerance_j: f64 },
}
impl From<JkrCombinedError> for JkrCombinedLoadingError {
    fn from(error: JkrCombinedError) -> Self {
        Self::Source(error)
    }
}
impl core::fmt::Display for JkrCombinedLoadingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "proportional combined JKR loading: {self:?}")
    }
}
impl std::error::Error for JkrCombinedLoadingError {}

/// Positive load magnitude at the first ideal boundary. The original material
/// envelope may refuse an earlier state; this is not a gross-slip threshold.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedLoadingLimit {
    pub kind: JkrCombinedBoundary,
    pub load_factor: f64,
    pub contact_radius_m: f64,
    /// Numerical root bracket, expanded for radius conversion roundoff.
    pub radius_bracket_m: [f64; 2],
    /// Original U_a at the rounded radius and actual rounded coordinates.
    pub radial_residual_n: f64,
    pub radial_residual_tolerance_n: f64,
    /// Original U_aa; near zero at a fold, generally positive at an opening edge.
    pub radial_curvature_n_per_m: f64,
    pub radial_curvature_tolerance_n_per_m: f64,
}

/// Accepted source-owned history. Public limit/response diagnostics cannot
/// manufacture an attached state or authorize its irreversible release.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedLoadingState {
    identity: ContentHash,
    source_identity: ContentHash,
    load_factor: f64,
    coordinates: [f64; 4],
    contact: Option<JkrCombinedResponse>,
    released_energy_j: f64,
    steps: u64,
}
impl JkrCombinedLoadingState {
    pub const fn model_identity(&self) -> ContentHash {
        self.source_identity
    }
    pub const fn load_factor(&self) -> f64 {
        self.load_factor
    }
    pub const fn coordinates(&self) -> [f64; 4] {
        self.coordinates
    }
    pub const fn contact(&self) -> Option<JkrCombinedResponse> {
        self.contact
    }
    pub const fn released_energy_j(&self) -> f64 {
        self.released_energy_j
    }
    pub const fn steps(&self) -> u64 {
        self.steps
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.contact.map_or(0.0, |point| point.reversible_energy_j)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedLoadingReport {
    pub before: Option<JkrCombinedResponse>,
    pub contact: Option<JkrCombinedResponse>,
    pub boundary: Option<JkrCombinedLoadingLimit>,
    pub requested_load_factor: f64,
    pub accepted_load_factor: f64,
    pub untraversed_load_factor: f64,
    pub normal_elastic_change_j: f64,
    pub tangential_elastic_change_j: f64,
    pub torsional_elastic_change_j: f64,
    pub surface_energy_change_j: f64,
    /// Integral T.du + M.dphi; the held normal coordinate does no work.
    pub prescribed_work_j: f64,
    /// Held-coordinate mechanical release; a downstream owner must route it.
    pub released_energy_j: f64,
    pub storage_change_j: f64,
    /// Delta U + release - prescribed work, independently checked at endpoints.
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub minimum_contact_radius_m: f64,
    pub maximum_contact_radius_m: f64,
    pub evaluations: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedLoadingStep {
    previous: JkrCombinedLoadingState,
    state: JkrCombinedLoadingState,
    report: JkrCombinedLoadingReport,
}
impl JkrCombinedLoadingStep {
    pub const fn previous_state(&self) -> JkrCombinedLoadingState {
        self.previous
    }
    pub const fn state(&self) -> JkrCombinedLoadingState {
        self.state
    }
    pub const fn report(&self) -> &JkrCombinedLoadingReport {
        &self.report
    }
}

#[derive(Clone, Copy, Debug)]
pub struct JkrCombinedLoading<'a> {
    source: &'a SourceJkrMixedMode,
    indentation_m: f64,
    direction: [f64; 3],
    controls: JkrCombinedLoadingControls,
    identity: ContentHash,
    // A detached thermal successor needs no hypothetical attached preload.
    zero: Option<JkrCombinedResponse>,
    limit: Option<JkrCombinedLoadingLimit>,
    preparation_evaluations: usize,
}

impl<'a> JkrCombinedLoading<'a> {
    /// Prepare a proportional loading family around an admitted zero-load
    /// contact. Pure twist at zero approach correctly has no finite boundary.
    /// An ideal limit outside the source envelope is retained diagnostically;
    /// reaching it still requires all original source admission checks.
    pub fn new(
        source: &'a SourceJkrMixedMode,
        indentation_m: f64,
        shear_direction_m: [f64; 2],
        twist_direction_rad: f64,
        controls: JkrCombinedLoadingControls,
        gate: &CancelGate,
    ) -> Result<Self, JkrCombinedLoadingError> {
        poll(gate)?;
        let direction = [
            shear_direction_m[0],
            shear_direction_m[1],
            twist_direction_rad,
        ];
        let mut loading = Self::unprepared(source, indentation_m, direction, controls)?;
        let mut work = Work::new(controls.maximum_evaluations, gate);
        loading.prepare(&mut work)?;
        poll(gate)?;
        Ok(loading)
    }

    fn unprepared(
        source: &'a SourceJkrMixedMode,
        indentation_m: f64,
        direction: [f64; 3],
        controls: JkrCombinedLoadingControls,
    ) -> Result<Self, JkrCombinedLoadingError> {
        if direction.iter().any(|x| !x.is_finite())
            || direction.iter().all(|x| *x == 0.0)
            || !direction[0].hypot(direction[1]).is_finite()
            || !indentation_m.is_finite()
            || controls.maximum_iterations == 0
            || controls.maximum_evaluations == 0
            || [
                controls.energy_absolute_tolerance_j,
                controls.energy_relative_tolerance,
            ]
            .into_iter()
            .any(|x| !x.is_finite() || x < 0.0)
        {
            return Err(JkrCombinedLoadingError::Invalid(
                "finite proportional loading controls",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-combined-loading.v1");
        hash.update(source.identity().as_bytes());
        for value in direction.into_iter().chain([
            indentation_m,
            controls.energy_absolute_tolerance_j,
            controls.energy_relative_tolerance,
        ]) {
            hash.update(&value.to_bits().to_le_bytes());
        }
        for value in [controls.maximum_iterations, controls.maximum_evaluations] {
            hash.update(&(value as u64).to_le_bytes());
        }
        Ok(Self {
            source,
            indentation_m,
            direction,
            controls,
            identity: hash.finalize(),
            zero: None,
            limit: None,
            preparation_evaluations: 0,
        })
    }

    fn prepare(&mut self, work: &mut Work<'_>) -> Result<(), JkrCombinedLoadingError> {
        let start = work.evaluations;
        work.charge()?;
        self.zero = Some(self.source.evaluate_combined_displacement_controlled([
            self.indentation_m,
            0.0,
            0.0,
            0.0,
        ])?);
        self.limit = self.locate_limit(work)?;
        self.preparation_evaluations = work.evaluations - start;
        Ok(())
    }

    pub const fn source(&self) -> &'a SourceJkrMixedMode {
        self.source
    }
    /// Prepared first boundary. `None` also applies to a detached thermal
    /// successor, which makes no attached-family claim and cannot initialize
    /// an attached state.
    pub const fn limit(&self) -> Option<JkrCombinedLoadingLimit> {
        self.limit
    }
    pub const fn controls(&self) -> JkrCombinedLoadingControls {
        self.controls
    }
    pub const fn preparation_evaluations(&self) -> usize {
        self.preparation_evaluations
    }

    /// Initialize an independently prescribed attached preload. An owner
    /// continued only as detached history has no such prepared family.
    pub fn initial_state(
        &self,
        load_factor: f64,
        gate: &CancelGate,
    ) -> Result<JkrCombinedLoadingState, JkrCombinedLoadingError> {
        let mut work = Work::new(self.controls.maximum_evaluations, gate);
        let contact = self.evaluate(load_factor, &mut work)?;
        poll(gate)?;
        Ok(JkrCombinedLoadingState {
            identity: self.identity,
            source_identity: self.source.identity(),
            load_factor,
            coordinates: contact.coordinates,
            contact: Some(contact),
            released_energy_j: 0.0,
            steps: 0,
        })
    }

    /// Traverse a signed straight loading segment up to its FIRST boundary.
    /// Elastic unloading and reversal return signed work without new release.
    /// Detached history remains detached, even when lambda returns to zero.
    pub fn advance_to_load_factor(
        &self,
        previous: JkrCombinedLoadingState,
        requested: f64,
        gate: &CancelGate,
    ) -> Result<JkrCombinedLoadingStep, JkrCombinedLoadingError> {
        self.admit(previous, gate)?;
        self.coordinates(requested)?;
        finite(requested - previous.load_factor, "finite loading increment")?;
        let mut work = Work::new(self.controls.maximum_evaluations, gate);
        let mut accepted = requested;
        let mut radii = [0.0; 2];
        let contact = if let Some(before) = previous.contact {
            if let Some(limit) = self.limit {
                accepted = requested.clamp(-limit.load_factor, limit.load_factor);
            }
            let after = if accepted == previous.load_factor {
                before
            } else {
                self.evaluate(accepted, &mut work)?
            };
            radii = [
                before.contact_radius_m.min(after.contact_radius_m),
                before.contact_radius_m.max(after.contact_radius_m),
            ];
            if previous.load_factor.is_sign_negative() != accepted.is_sign_negative() {
                radii[1] = radii[1].max(self.prepared_zero()?.contact_radius_m);
            }
            radii = [radii[0].next_down(), radii[1].next_up()];
            for a in radii {
                self.source
                    .normal
                    .check_applicability(a)
                    .map_err(JkrCombinedError::from)?;
            }
            let maximum = previous.load_factor.abs().max(accepted.abs());
            self.source.admit_combined_nominal_bounds(
                maximum * self.direction[0].hypot(self.direction[1]),
                maximum * self.direction[2].abs(),
                radii[0],
            )?;
            Some(after)
        } else {
            None
        };
        let state = JkrCombinedLoadingState {
            load_factor: accepted,
            coordinates: self.coordinates(accepted)?,
            contact,
            steps: next_step(previous.steps)?,
            ..previous
        };
        let changes = self.storage_changes(previous.contact, contact)?;
        let prescribed = sum(changes)?;
        self.finish(
            previous,
            state,
            requested,
            changes,
            prescribed,
            0.0,
            radii,
            work.evaluations,
            gate,
        )
    }

    /// Release the complete original potential only at a resolved radial fold
    /// with negative held approach. There is no admissible separated landing
    /// at an opening edge or nonnegative approach. Repeated release refuses.
    pub fn detach_at_fold(
        &self,
        previous: JkrCombinedLoadingState,
        gate: &CancelGate,
    ) -> Result<JkrCombinedLoadingStep, JkrCombinedLoadingError> {
        self.admit(previous, gate)?;
        let before = previous
            .contact
            .filter(|point| {
                point.boundary == Some(JkrCombinedBoundary::RadialFold)
                    && self.indentation_m < 0.0
                    && self
                        .limit
                        .is_some_and(|limit| previous.load_factor.abs() == limit.load_factor)
            })
            .ok_or(JkrCombinedLoadingError::NotAtSeparatingFold)?;
        let release = finite(before.reversible_energy_j, "finite fold release")?;
        if release < 0.0 {
            return Err(JkrCombinedLoadingError::Invalid("nonnegative fold release"));
        }
        let total = finite(
            (dd(previous.released_energy_j) + dd(release)).to_f64(),
            "finite cumulative release",
        )?;
        let state = JkrCombinedLoadingState {
            contact: None,
            released_energy_j: total,
            steps: next_step(previous.steps)?,
            ..previous
        };
        let changes = self.storage_changes(previous.contact, None)?;
        self.finish(
            previous,
            state,
            previous.load_factor,
            changes,
            0.0,
            release,
            [0.0, before.contact_radius_m],
            0,
            gate,
        )
    }

    fn admit(
        &self,
        state: JkrCombinedLoadingState,
        gate: &CancelGate,
    ) -> Result<(), JkrCombinedLoadingError> {
        poll(gate)?;
        if state.identity != self.identity {
            return Err(JkrCombinedLoadingError::StateMismatch);
        }
        Ok(())
    }

    fn coordinates(&self, load: f64) -> Result<[f64; 4], JkrCombinedLoadingError> {
        let q = [
            self.indentation_m,
            load * self.direction[0],
            load * self.direction[1],
            load * self.direction[2],
        ];
        if !load.is_finite() || q.iter().any(|x| !x.is_finite()) {
            return Err(JkrCombinedLoadingError::Invalid(
                "finite prescribed coordinates",
            ));
        }
        Ok(q)
    }

    fn evaluate(
        &self,
        load: f64,
        work: &mut Work<'_>,
    ) -> Result<JkrCombinedResponse, JkrCombinedLoadingError> {
        poll(work.gate)?;
        let zero = self.prepared_zero()?;
        let q = self.coordinates(load)?;
        if let Some(limit) = self.limit {
            if load.abs() > limit.load_factor {
                return Err(JkrCombinedLoadingError::BeyondFirstBoundary);
            }
            if load.abs() == limit.load_factor {
                work.charge()?;
                let curvature = if limit.kind == JkrCombinedBoundary::RadialFold {
                    0.0
                } else {
                    limit.radial_curvature_n_per_m.max(0.0)
                };
                return Ok(self.source.combined_response_at_radius(
                    q,
                    limit.contact_radius_m,
                    curvature,
                    Some(limit.kind),
                )?);
            }
        }
        if load == 0.0 {
            return Ok(zero);
        }
        work.charge()?;
        Ok(self.source.evaluate_combined_displacement_controlled(q)?)
    }

    fn locate_limit(
        &self,
        work: &mut Work<'_>,
    ) -> Result<Option<JkrCombinedLoadingLimit>, JkrCombinedLoadingError> {
        let zero = self.prepared_zero()?;
        if zero.boundary == Some(JkrCombinedBoundary::RadialFold) {
            return self
                .make_limit(
                    JkrCombinedBoundary::RadialFold,
                    0.0,
                    zero.contact_radius_m,
                    [zero.contact_radius_m; 2],
                )
                .map(Some);
        }
        let normal = self.source.normal_source();
        let scale = dd(normal.spinodal_radius_m);
        let depth = scale * scale / dd(normal.effective_radius_m());
        let e = dd(normal.normal_state().elastic().reduced_modulus_pa());
        let stationary_scale = e * depth * depth;
        let d = dd(self.indentation_m) / depth;
        let b = dd(8.0)
            * dd(self.source.torsional_reduced_shear_modulus_pa())
            * scale
            * scale
            * dd(self.direction[2])
            * dd(self.direction[2])
            / stationary_scale;
        let shear_square = dd(self.direction[0]) * dd(self.direction[0])
            + dd(self.direction[1]) * dd(self.direction[1]);
        let t = dd(0.5) * dd(self.source.tangential_stiffness_per_radius) * shear_square
            / stationary_scale;
        let ell = dd(2.0) * dd(PI) * dd(normal.normal_state().adhesion_energy_j_per_m2()) * scale
            / stationary_scale;
        let boundary = radius::locate(
            d,
            b,
            t,
            ell,
            zero.contact_radius_m / normal.spinodal_radius_m,
            self.controls.maximum_iterations,
            || work.charge(),
        )
        .map_err(|error| match error {
            radius::Error::Probe(error) => error,
            radius::Error::IterationBudget => JkrCombinedLoadingError::Budget,
            radius::Error::Unresolved => JkrCombinedLoadingError::UnresolvedBoundary,
        })?;
        let Some(boundary) = boundary else {
            return Ok(None);
        };
        work.charge()?;
        let mut a = normal.spinodal_radius_m * boundary.radius;
        let kind = match boundary.kind {
            radius::Kind::Fold => JkrCombinedBoundary::RadialFold,
            radius::Kind::Opening => JkrCombinedBoundary::OpeningEdge,
        };
        // Keep the rounded opening-edge diagnostic on the admitted side.
        if kind == JkrCombinedBoundary::OpeningEdge
            && a * (a / normal.effective_radius_m()) < self.indentation_m
        {
            a = a.next_up();
        }
        let y = dd(a) / scale;
        let gap = y * y - d;
        let square = (ell * y - gap * gap) / (b * y * y + t);
        if !square.hi.is_finite() || !square.lo.is_finite() || square.lt(Dd::ZERO) {
            return Err(JkrCombinedLoadingError::UnresolvedBoundary);
        }
        let load = square.to_f64().sqrt();
        let bracket = [
            (normal.spinodal_radius_m * boundary.bracket[0]).next_down(),
            (normal.spinodal_radius_m * boundary.bracket[1])
                .max(a)
                .next_up(),
        ];
        self.make_limit(kind, load, a, bracket).map(Some)
    }

    fn prepared_zero(&self) -> Result<JkrCombinedResponse, JkrCombinedLoadingError> {
        self.zero.ok_or(JkrCombinedLoadingError::Invalid(
            "detached continuation has no prepared attached loading family",
        ))
    }

    fn make_limit(
        &self,
        kind: JkrCombinedBoundary,
        load: f64,
        a: f64,
        bracket: [f64; 2],
    ) -> Result<JkrCombinedLoadingLimit, JkrCombinedLoadingError> {
        let q = self.coordinates(load)?;
        let normal = self.source.normal_source();
        let e = normal.normal_state().elastic().reduced_modulus_pa();
        let c = 2.0 * PI * normal.normal_state().adhesion_energy_j_per_m2();
        let opening = a * (a / normal.effective_radius_m()) - q[0];
        let shear = 0.5 * self.source.tangential_stiffness_per_radius * (q[1] * q[1] + q[2] * q[2]);
        let twist = 8.0 * self.source.torsional_reduced_shear_modulus_pa() * q[3] * q[3];
        let radial_scale = e * opening * opening + shear + twist * a * a + c * a;
        let curvature_scale =
            4.0 * e * a / normal.effective_radius_m() * opening + 2.0 * twist * a + c;
        let residual = self.source.combined_radial_residual(q, a);
        let curvature = self.source.combined_radial_curvature(q, a);
        let radial_tolerance = 256.0 * f64::EPSILON * radial_scale;
        let curvature_tolerance = 256.0 * f64::EPSILON * curvature_scale;
        if [
            load,
            a,
            opening,
            residual,
            curvature,
            radial_tolerance,
            curvature_tolerance,
        ]
        .into_iter()
        .any(|x| !x.is_finite())
            || a <= 0.0
            || load < 0.0
            || opening < 0.0
            || residual.abs() > radial_tolerance
            || (kind == JkrCombinedBoundary::RadialFold && curvature.abs() > curvature_tolerance)
            || (kind == JkrCombinedBoundary::OpeningEdge && curvature < -curvature_tolerance)
        {
            return Err(JkrCombinedLoadingError::UnresolvedBoundary);
        }
        Ok(JkrCombinedLoadingLimit {
            kind,
            load_factor: load,
            contact_radius_m: a,
            radius_bracket_m: bracket,
            radial_residual_n: residual,
            radial_residual_tolerance_n: radial_tolerance,
            radial_curvature_n_per_m: curvature,
            radial_curvature_tolerance_n_per_m: curvature_tolerance,
        })
    }

    // Factor the ORIGINAL polynomial potential; no singular equilibrium
    // tangent, numerical quadrature or subtraction of large total energies.
    fn storage_changes(
        &self,
        before: Option<JkrCombinedResponse>,
        after: Option<JkrCombinedResponse>,
    ) -> Result<[f64; 4], JkrCombinedLoadingError> {
        let (Some(before), Some(after)) = (before, after) else {
            return match (before, after) {
                (None, None) => Ok([0.0; 4]),
                (Some(point), None) => Ok([
                    -point.normal_elastic_energy_j,
                    -point.tangential_elastic_energy_j,
                    -point.torsional_elastic_energy_j,
                    -point.surface_energy_j,
                ]),
                _ => Err(JkrCombinedLoadingError::Invalid("no unowned rebonding")),
            };
        };
        if before == after {
            return Ok([0.0; 4]);
        }
        let a0 = dd(before.contact_radius_m);
        let a1 = dd(after.contact_radius_m);
        let da = a1 - a0;
        let a0sq = a0 * a0;
        let a1sq = a1 * a1;
        let cube_difference = da * (a1sq + a1 * a0 + a0sq);
        let fifth_difference =
            da * (a1sq * a1sq + a1sq * a1 * a0 + a1sq * a0sq + a1 * a0sq * a0 + a0sq * a0sq);
        let normal = self.source.normal_source();
        let delta = dd(self.indentation_m);
        let r = dd(normal.effective_radius_m());
        let normal_change = dd(normal.normal_state().elastic().reduced_modulus_pa())
            * (da * delta * delta - dd(2.0) * delta * cube_difference / (dd(3.0) * r)
                + fifth_difference / (dd(5.0) * r * r));
        let mut shear_change = Dd::ZERO;
        for i in 1..3 {
            let u0 = dd(before.coordinates[i]);
            let u1 = dd(after.coordinates[i]);
            shear_change = shear_change
                + da * (u1 * u1 + u0 * u0) / dd(2.0)
                + (a1 + a0) / dd(2.0) * (u1 - u0) * (u1 + u0);
        }
        shear_change = shear_change * dd(0.5) * dd(self.source.tangential_stiffness_per_radius);
        let p0 = dd(before.coordinates[3]);
        let p1 = dd(after.coordinates[3]);
        let torsional_change = dd(8.0) / dd(3.0)
            * dd(self.source.torsional_reduced_shear_modulus_pa())
            * (cube_difference * (p1 * p1 + p0 * p0) / dd(2.0)
                + (a1 * a1sq + a0 * a0sq) / dd(2.0) * (p1 - p0) * (p1 + p0));
        let surface_change =
            -dd(PI) * dd(normal.normal_state().adhesion_energy_j_per_m2()) * da * (a1 + a0);
        let changes = [
            normal_change,
            shear_change,
            torsional_change,
            surface_change,
        ]
        .map(Dd::to_f64);
        for value in changes {
            finite(value, "finite original storage changes")?;
        }
        Ok(changes)
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        previous: JkrCombinedLoadingState,
        state: JkrCombinedLoadingState,
        requested: f64,
        changes: [f64; 4],
        prescribed: f64,
        release: f64,
        radii: [f64; 2],
        evaluations: usize,
        gate: &CancelGate,
    ) -> Result<JkrCombinedLoadingStep, JkrCombinedLoadingError> {
        let change = finite(
            (dd(state.stored_energy_j()) - dd(previous.stored_energy_j())).to_f64(),
            "finite storage change",
        )?;
        let residual = finite(
            (dd(change) + dd(release) - dd(prescribed)).to_f64(),
            "finite loading balance",
        )?;
        let mut scale = prescribed.abs() + release;
        for point in [previous.contact, state.contact].into_iter().flatten() {
            scale += point.normal_elastic_energy_j
                + point.tangential_elastic_energy_j
                + point.torsional_elastic_energy_j
                - point.surface_energy_j;
        }
        let tolerance = finite(
            self.controls.energy_absolute_tolerance_j
                + self.controls.energy_relative_tolerance * scale,
            "finite loading tolerance",
        )?;
        if residual.abs() > tolerance {
            return Err(JkrCombinedLoadingError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            });
        }
        poll(gate)?;
        Ok(JkrCombinedLoadingStep {
            previous,
            state,
            report: JkrCombinedLoadingReport {
                before: previous.contact,
                contact: state.contact,
                boundary: if state.contact.is_some_and(|point| point.boundary.is_some())
                    || release != 0.0
                {
                    self.limit
                } else {
                    None
                },
                requested_load_factor: requested,
                accepted_load_factor: state.load_factor,
                untraversed_load_factor: requested - state.load_factor,
                normal_elastic_change_j: changes[0],
                tangential_elastic_change_j: changes[1],
                torsional_elastic_change_j: changes[2],
                surface_energy_change_j: changes[3],
                prescribed_work_j: prescribed,
                released_energy_j: release,
                storage_change_j: change,
                energy_residual_j: residual,
                energy_tolerance_j: tolerance,
                minimum_contact_radius_m: radii[0],
                maximum_contact_radius_m: radii[1],
                evaluations,
            },
        })
    }
}

struct Work<'a> {
    gate: &'a CancelGate,
    maximum: usize,
    evaluations: usize,
}
impl<'a> Work<'a> {
    fn new(maximum: usize, gate: &'a CancelGate) -> Self {
        Self {
            gate,
            maximum,
            evaluations: 0,
        }
    }
    fn charge(&mut self) -> Result<(), JkrCombinedLoadingError> {
        poll(self.gate)?;
        if self.evaluations >= self.maximum {
            return Err(JkrCombinedLoadingError::Budget);
        }
        self.evaluations += 1;
        Ok(())
    }
}
fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}
fn sum(values: [f64; 4]) -> Result<f64, JkrCombinedLoadingError> {
    finite(
        values
            .into_iter()
            .fold(Dd::ZERO, |total, value| total + dd(value))
            .to_f64(),
        "finite prescribed work",
    )
}
fn finite(value: f64, field: &'static str) -> Result<f64, JkrCombinedLoadingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(JkrCombinedLoadingError::Invalid(field))
    }
}
fn next_step(steps: u64) -> Result<u64, JkrCombinedLoadingError> {
    steps
        .checked_add(1)
        .ok_or(JkrCombinedLoadingError::Invalid("loading step overflow"))
}
fn poll(gate: &CancelGate) -> Result<(), JkrCombinedLoadingError> {
    if gate.is_requested() {
        Err(JkrCombinedLoadingError::Cancelled)
    } else {
        Ok(())
    }
}
