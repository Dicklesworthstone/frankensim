//! Conservative nonlinear evolution of a periodic Newtonian slender jet.
//!
//! This advances the axisymmetric one-dimensional equations of Eggers and
//! Dupont (1994), equations (17)-(18), with full surface-of-revolution curvature:
//! <https://people.maths.bris.ac.uk/~majge/ED94.pdf>.
//!
//! With A=h^2 and q=A*u, the equations without gravity are
//! A_t=-q_z and q_t=-[q*u+(sigma/rho)*(A*kappa-2*h/sqrt(1+h_z^2))
//! -3*nu*A*u_z]_z, where kappa=1/(h*sqrt(1+h_z^2))
//! -h_zz/(1+h_z^2)^(3/2). The capillary stress uses the identity
//! A*kappa_z=(A*kappa-2*h/sqrt(1+h_z^2))_z.
//!
//! A is a cell value; q and u live at faces i+1/2. Centered curvature and
//! viscosity preserve the circular geometry. Mass-consistent upwind
//! momentum transport uses a minmod-limited linear reconstruction and exposes
//! its numerical dissipation. Explicit SSPRK3 is the default; optional IMEX
//! treats viscosity implicitly. Adaptive steps retain advection, drainage and
//! capillary-Jacobian bounds; explicit mode also limits viscous diffusion. Smooth transport, capillary and viscous
//! terms are second order in space; the limiter reduces order at extrema.
//!
//! All coefficients are constant. The liquid is incompressible/Newtonian,
//! its exterior has negligible inertia/viscosity, and the domain is periodic.
//! No nozzle, stretching imposed by gravity, heat, gas drag, breakup topology,
//! coalescence, viscoelasticity or post-pinch continuation is represented.
//! Authored radius, slope and resolution limits stop the calculation before
//! leaving the admitted model. They are not error certificates.
//!
//! Eleven reusable f64 arrays give O(N) storage in explicit mode. Selecting
//! implicit viscosity prepares twelve additional arrays once; stepping allocates nothing.
//! Fixed serial order and arithmetic/sqrt give deterministic evaluations.
//! Cancellation is polled every 256 cells and immediately before publication.

use fs_exec::CancelGate;

mod implicit;
pub use implicit::{SlenderJetIntegrationMethod, SlenderJetTemporalError};

/// Constant SI properties of an explicitly Newtonian liquid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlenderJetParameters {
    pub density_kg_m3: f64,
    pub surface_tension_n_m: f64,
    pub kinematic_viscosity_m2_s: f64,
}

/// Author-selected boundaries of the resolved slender-jet model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlenderJetLimits {
    /// Positive neck radius below which this calculation must stop [m].
    pub minimum_radius_m: f64,
    /// Maximum absolute adjacent-cell radius slope; must lie in (0, 1].
    pub maximum_abs_slope: f64,
    /// Minimum h/dz, at least two; no radius floor or implicit remeshing occurs.
    pub minimum_radius_cells: f64,
}

/// Explicit accuracy/work controls for one advance call.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlenderJetStepBudget {
    /// Positive upper timestep in seconds, independent of stability limits.
    pub maximum_step_s: f64,
    /// Safety factor on the combined rate bound, in (0, 0.25].
    pub cfl: f64,
    /// Maximum attempted substeps, including retries after a changing rate bound.
    pub maximum_substeps: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SlenderJetLimitKind {
    MinimumRadius,
    MaximumSlope,
    RadiusResolution,
}

/// An authored physical/resolution limit exceeded by an uncommitted stage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlenderJetViolation {
    pub kind: SlenderJetLimitKind,
    pub cell: usize,
    /// Radius [m], absolute slope, or h/dz, according to kind.
    pub observed: f64,
    pub limit: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SlenderJetOutcome {
    Completed,
    /// The last admitted substep is published. This is not an exact event time
    /// or a claim that pinch-off has occurred; the rejected RK stage is absent.
    DomainLimit(SlenderJetViolation),
}

/// Diagnostics of accepted state, with physical factors pi and rho included.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlenderJetDiagnostics {
    pub volume_m3: f64,
    pub axial_momentum_kg_m_s: f64,
    /// Momentum divided by liquid mass; no reference-frame shift is applied.
    pub mean_velocity_m_s: f64,
    pub minimum_radius_m: f64,
    pub maximum_radius_m: f64,
    pub maximum_abs_slope: f64,
    pub kinetic_energy_j: f64,
    /// Area of the piecewise conical surface times sigma.
    pub surface_energy_j: f64,
}

/// Advance result. Numerical failure/cancellation returns Err instead and
/// leaves the entire pre-call state, diagnostics and time unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlenderJetAdvance {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub requested_duration_s: f64,
    pub accepted_substeps: usize,
    pub attempted_substeps: usize,
    pub outcome: SlenderJetOutcome,
    pub before: SlenderJetDiagnostics,
    pub after: SlenderJetDiagnostics,
    /// Numerical integration policy used for this advance.
    pub integration_method: SlenderJetIntegrationMethod,
    /// Largest embedded local time-error indicator across accepted substeps.
    /// Present for implicit viscosity; not a global error bound or certificate.
    pub temporal_error: Option<SlenderJetTemporalError>,
    /// RK quadrature of 3*pi*rho*nu*integral(A*u_z^2 dz), not evolved heat.
    pub viscous_dissipation_j: f64,
    /// RK quadrature of limited-upwind momentum dissipation; not material viscosity.
    pub upwind_dissipation_j: f64,
    /// E_after-E_before plus both accounted losses. Signed truncation/roundoff
    /// remainder; the full nonlinear discretization is not energy exact.
    pub energy_balance_residual_j: f64,
}

/// Admission/numerical refusal; advance never partially publishes an Err.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NonlinearJetError {
    InvalidInput(&'static str),
    InvalidState(&'static str),
    NonRepresentable,
    Allocation,
    Cancelled,
    SubstepBudget,
    /// Constructor refusal. During advance this becomes a normal DomainLimit.
    InitialDomain(SlenderJetViolation),
}
impl core::fmt::Display for NonlinearJetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "periodic slender jet: {self:?}")
    }
}
impl std::error::Error for NonlinearJetError {}

/// Shared evolution interface for native and material-bound jet owners.
pub trait SlenderJetEvolution {
    /// Return the accepted state whose constant coefficients this owner retains.
    fn state(&self) -> &PeriodicSlenderJet;

    /// Advance that same state. Errors preserve the full pre-call state; a
    /// domain-limit report retains only completed, admitted substeps. The
    /// returned report and state must describe the same physical trajectory.
    fn advance(
        &mut self,
        duration_s: f64,
        budget: SlenderJetStepBudget,
        gate: &CancelGate,
    ) -> Result<SlenderJetAdvance, NonlinearJetError>;
}

#[derive(Clone, Debug)]
struct Fields {
    a: Vec<f64>,
    q: Vec<f64>,
}
impl Fields {
    fn zeros(n: usize, gate: &CancelGate) -> Result<Self, NonlinearJetError> {
        Ok(Self {
            a: zeros(n, gate)?,
            q: zeros(n, gate)?,
        })
    }
}
#[derive(Clone, Debug)]
struct Spatial {
    h: Vec<f64>,
    u: Vec<f64>,
    flux: Vec<f64>,
}
#[derive(Clone, Copy)]
struct Evaluation {
    diagnostics: SlenderJetDiagnostics,
    rate_per_s: f64,
    viscous_power_w: f64,
    upwind_power_w: f64,
}

/// Uniform periodic grid with a transaction workspace, all in SI units.
#[derive(Clone, Debug)]
pub struct PeriodicSlenderJet {
    parameters: SlenderJetParameters,
    limits: SlenderJetLimits,
    length_m: f64,
    spacing_m: f64,
    time_s: f64,
    accepted: Fields,
    candidate: Fields,
    stage: Fields,
    rhs: Fields,
    spatial: Spatial,
    diagnostics: SlenderJetDiagnostics,
    integration_method: SlenderJetIntegrationMethod,
    implicit_workspace: Option<implicit::Workspace>,
}

impl PeriodicSlenderJet {
    /// Construct at t=0 from N cell radii and N FACE velocities at i+1/2.
    /// N>=8. The axis is periodic with cell centers z=i*length/N. Physical
    /// cross-sectional area is pi*h^2 and physical volume flux is pi*q.
    pub fn new(
        parameters: SlenderJetParameters,
        length_m: f64,
        radii_m: &[f64],
        velocities_m_s: &[f64],
        limits: SlenderJetLimits,
        gate: &CancelGate,
    ) -> Result<Self, NonlinearJetError> {
        poll(gate)?;
        validate_parameters(parameters, limits)?;
        let n = radii_m.len();
        if n < 8 || velocities_m_s.len() != n {
            return Err(NonlinearJetError::InvalidInput("grid/profile lengths"));
        }
        if !length_m.is_finite() || length_m <= 0.0 {
            return Err(NonlinearJetError::InvalidInput("period length"));
        }
        let dx = positive(length_m / n as f64)?;
        let mut accepted = Fields::zeros(n, gate)?;
        for i in 0..n {
            poll_at(gate, i)?;
            let h = radii_m[i];
            if !h.is_finite() || h <= 0.0 || !velocities_m_s[i].is_finite() {
                return Err(NonlinearJetError::InvalidInput("initial profile"));
            }
            accepted.a[i] = positive(h * h)?;
        }
        for (i, &u) in velocities_m_s.iter().enumerate() {
            poll_at(gate, i)?;
            let q = average(accepted.a[i], accepted.a[next(i, n)]) * u;
            if !q.is_finite() || (u != 0.0 && q == 0.0) {
                return Err(NonlinearJetError::NonRepresentable);
            }
            accepted.q[i] = q;
        }
        let candidate = Fields::zeros(n, gate)?;
        let stage = Fields::zeros(n, gate)?;
        let mut rhs = Fields::zeros(n, gate)?;
        let mut spatial = Spatial {
            h: zeros(n, gate)?,
            u: zeros(n, gate)?,
            flux: zeros(n, gate)?,
        };
        let initial = evaluate(
            &accepted,
            parameters,
            limits,
            dx,
            &mut spatial,
            &mut rhs,
            gate,
        )?;
        Ok(Self {
            parameters,
            limits,
            length_m,
            spacing_m: dx,
            time_s: 0.0,
            accepted,
            candidate,
            stage,
            rhs,
            spatial,
            diagnostics: initial.diagnostics,
            integration_method: SlenderJetIntegrationMethod::ExplicitSspRk3,
            implicit_workspace: None,
        })
    }

    #[must_use]
    pub const fn parameters(&self) -> SlenderJetParameters {
        self.parameters
    }
    #[must_use]
    pub const fn limits(&self) -> SlenderJetLimits {
        self.limits
    }
    #[must_use]
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    #[must_use]
    pub const fn spacing_m(&self) -> f64 {
        self.spacing_m
    }
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Cell h^2 values [m^2], NOT physical cross-sectional areas.
    #[must_use]
    pub fn radius_squared_m2(&self) -> &[f64] {
        &self.accepted.a
    }
    /// Face q=(h^2)_face*u [m^3/s], NOT physical volume fluxes pi*q.
    #[must_use]
    pub fn radius_squared_flux_m3_s(&self) -> &[f64] {
        &self.accepted.q
    }
    pub fn radii_m(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.accepted.a.iter().map(|a| a.sqrt())
    }
    /// Face velocities at i+1/2, in the same frame supplied at construction.
    pub fn velocities_m_s(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        let n = self.accepted.a.len();
        self.accepted
            .q
            .iter()
            .enumerate()
            .map(move |(i, q)| q / average(self.accepted.a[i], self.accepted.a[next(i, n)]))
    }
    #[must_use]
    pub const fn diagnostics(&self) -> SlenderJetDiagnostics {
        self.diagnostics
    }

    /// Advance toward a requested duration using the selected integration policy.
    ///
    /// The default is explicit SSPRK3. ImplicitViscosity uses third-order
    /// IMEX-SSP3(4,3,3), retaining capillary/advection limits and an embedded
    /// local time-error indicator; maximum_step_s still controls accuracy.
    ///
    /// A domain stop publishes only completed, admitted substeps. Its last
    /// accepted time brackets neither an exact pinch event nor a precise
    /// continuum limit crossing. Errors (including work-budget exhaustion and
    /// cancellation) roll the entire call back; workspace contents may change.
    pub fn advance(
        &mut self,
        duration_s: f64,
        budget: SlenderJetStepBudget,
        gate: &CancelGate,
    ) -> Result<SlenderJetAdvance, NonlinearJetError> {
        poll(gate)?;
        if !duration_s.is_finite() || duration_s < 0.0 {
            return Err(NonlinearJetError::InvalidInput("duration"));
        }
        if !budget.maximum_step_s.is_finite()
            || budget.maximum_step_s <= 0.0
            || !budget.cfl.is_finite()
            || budget.cfl <= 0.0
            || budget.cfl > 0.25
            || budget.maximum_substeps == 0
        {
            return Err(NonlinearJetError::InvalidInput("step budget"));
        }
        let target_time = self.time_s + duration_s;
        if !target_time.is_finite() || (duration_s > 0.0 && target_time <= self.time_s) {
            return Err(NonlinearJetError::NonRepresentable);
        }
        let before = self.diagnostics;
        let explicit_viscosity =
            self.integration_method == SlenderJetIntegrationMethod::ExplicitSspRk3;
        let mut temporal_error = (!explicit_viscosity).then(SlenderJetTemporalError::default);
        let n = self.accepted.a.len();
        for i in 0..n {
            poll_at(gate, i)?;
            self.candidate.a[i] = self.accepted.a[i];
            self.candidate.q[i] = self.accepted.q[i];
        }
        let mut base = evaluate_split(
            &self.candidate,
            self.parameters,
            self.limits,
            self.spacing_m,
            &mut self.spatial,
            &mut self.rhs,
            explicit_viscosity,
            gate,
        )?;
        let (mut elapsed, mut step_cap) = (0.0, budget.maximum_step_s);
        let (mut accepted_substeps, mut attempted_substeps) = (0, 0);
        let (mut viscous_loss, mut upwind_loss) = (0.0, 0.0);
        let mut outcome = SlenderJetOutcome::Completed;
        while elapsed < duration_s {
            poll(gate)?;
            if attempted_substeps >= budget.maximum_substeps {
                return Err(NonlinearJetError::SubstepBudget);
            }
            let remaining = duration_s - elapsed;
            let dt = remaining
                .min(step_cap)
                .min(stable_step(base.rate_per_s, budget.cfl));
            positive(dt)?;
            let next_elapsed = if dt == remaining {
                duration_s
            } else {
                elapsed + dt
            };
            if next_elapsed <= elapsed || self.time_s + next_elapsed <= self.time_s + elapsed {
                return Err(NonlinearJetError::NonRepresentable);
            }
            attempted_substeps += 1;
            let result = if explicit_viscosity {
                for i in 0..n {
                    poll_at(gate, i)?;
                    self.stage.a[i] = self.candidate.a[i] + dt * self.rhs.a[i];
                    self.stage.q[i] = self.candidate.q[i] + dt * self.rhs.q[i];
                }
                self.finish_stages(dt, budget.cfl, base, gate)
                    .map(|result| {
                        result.map(|(end, physical, numerical)| (end, physical, numerical, None))
                    })
            } else {
                implicit::step(
                    &self.candidate,
                    &mut self.stage,
                    &mut self.rhs,
                    &mut self.spatial,
                    self.implicit_workspace
                        .as_mut()
                        .ok_or(NonlinearJetError::InvalidState(
                            "implicit workspace is absent",
                        ))?,
                    self.parameters,
                    self.limits,
                    self.spacing_m,
                    dt,
                    budget.cfl,
                    gate,
                )
            };
            let (end, physical_loss, numerical_loss, step_error) = match result {
                Err(NonlinearJetError::InitialDomain(violation)) => {
                    outcome = SlenderJetOutcome::DomainLimit(violation);
                    break;
                }
                Err(error) => return Err(error),
                Ok(None) => {
                    step_cap = positive(0.5 * dt)?;
                    // Rebuild the current-state RHS overwritten by the rejected stage.
                    base = evaluate_split(
                        &self.candidate,
                        self.parameters,
                        self.limits,
                        self.spacing_m,
                        &mut self.spatial,
                        &mut self.rhs,
                        explicit_viscosity,
                        gate,
                    )?;
                    continue;
                }
                Ok(Some(result)) => result,
            };
            if let (Some(maximum), Some(estimate)) = (&mut temporal_error, step_error) {
                maximum.maximum_radius_error_m = maximum
                    .maximum_radius_error_m
                    .max(estimate.maximum_radius_error_m);
                maximum.maximum_velocity_error_m_s = maximum
                    .maximum_velocity_error_m_s
                    .max(estimate.maximum_velocity_error_m_s);
            }
            viscous_loss = finite(viscous_loss + physical_loss)?;
            upwind_loss = finite(upwind_loss + numerical_loss)?;
            core::mem::swap(&mut self.candidate, &mut self.stage);
            base = end; // final evaluation already left the next RHS in the workspace
            elapsed = next_elapsed;
            accepted_substeps += 1;
            step_cap = budget.maximum_step_s;
        }
        let after = base.diagnostics;
        let residual = finite(
            (after.kinetic_energy_j - before.kinetic_energy_j)
                + (after.surface_energy_j - before.surface_energy_j)
                + viscous_loss
                + upwind_loss,
        )?;
        let end_time_s = if elapsed == duration_s {
            target_time
        } else {
            self.time_s + elapsed
        };
        let report = SlenderJetAdvance {
            start_time_s: self.time_s,
            end_time_s,
            requested_duration_s: duration_s,
            accepted_substeps,
            attempted_substeps,
            outcome,
            before,
            after,
            integration_method: self.integration_method,
            temporal_error,
            viscous_dissipation_j: viscous_loss,
            upwind_dissipation_j: upwind_loss,
            energy_balance_residual_j: residual,
        };
        poll(gate)?;
        core::mem::swap(&mut self.accepted, &mut self.candidate);
        self.time_s = end_time_s;
        self.diagnostics = after;
        Ok(report)
    }

    // None means the evolving rate bound requires a smaller step.
    fn finish_stages(
        &mut self,
        dt: f64,
        cfl: f64,
        base: Evaluation,
        gate: &CancelGate,
    ) -> Result<Option<(Evaluation, f64, f64)>, NonlinearJetError> {
        let first = evaluate(
            &self.stage,
            self.parameters,
            self.limits,
            self.spacing_m,
            &mut self.spatial,
            &mut self.rhs,
            gate,
        )?;
        if dt > stable_step(first.rate_per_s, cfl) {
            return Ok(None);
        }
        blend(&mut self.stage, &self.candidate, &self.rhs, dt, 0.25, gate)?;
        let second = evaluate(
            &self.stage,
            self.parameters,
            self.limits,
            self.spacing_m,
            &mut self.spatial,
            &mut self.rhs,
            gate,
        )?;
        if dt > stable_step(second.rate_per_s, cfl) {
            return Ok(None);
        }
        blend(
            &mut self.stage,
            &self.candidate,
            &self.rhs,
            dt,
            2.0 / 3.0,
            gate,
        )?;
        let end = evaluate(
            &self.stage,
            self.parameters,
            self.limits,
            self.spacing_m,
            &mut self.spatial,
            &mut self.rhs,
            gate,
        )?;
        if dt > stable_step(end.rate_per_s, cfl) {
            return Ok(None);
        }
        let viscous = finite(
            dt * (base.viscous_power_w / 6.0
                + first.viscous_power_w / 6.0
                + (2.0 / 3.0) * second.viscous_power_w),
        )?;
        let upwind = finite(
            dt * (base.upwind_power_w / 6.0
                + first.upwind_power_w / 6.0
                + (2.0 / 3.0) * second.upwind_power_w),
        )?;
        Ok(Some((end, viscous, upwind)))
    }
}

impl SlenderJetEvolution for PeriodicSlenderJet {
    fn state(&self) -> &PeriodicSlenderJet {
        self
    }
    fn advance(
        &mut self,
        duration_s: f64,
        budget: SlenderJetStepBudget,
        gate: &CancelGate,
    ) -> Result<SlenderJetAdvance, NonlinearJetError> {
        PeriodicSlenderJet::advance(self, duration_s, budget, gate)
    }
}

fn validate_parameters(
    p: SlenderJetParameters,
    limits: SlenderJetLimits,
) -> Result<(), NonlinearJetError> {
    if !p.density_kg_m3.is_finite()
        || p.density_kg_m3 <= 0.0
        || !p.surface_tension_n_m.is_finite()
        || p.surface_tension_n_m < 0.0
        || !p.kinematic_viscosity_m2_s.is_finite()
        || p.kinematic_viscosity_m2_s < 0.0
    {
        return Err(NonlinearJetError::InvalidInput("liquid coefficients"));
    }
    if !limits.minimum_radius_m.is_finite()
        || limits.minimum_radius_m <= 0.0
        || !limits.maximum_abs_slope.is_finite()
        || limits.maximum_abs_slope <= 0.0
        || limits.maximum_abs_slope > 1.0
        || !limits.minimum_radius_cells.is_finite()
        || limits.minimum_radius_cells < 2.0
    {
        return Err(NonlinearJetError::InvalidInput("slender/resolution limits"));
    }
    if p.surface_tension_n_m > 0.0 {
        positive(p.surface_tension_n_m / p.density_kg_m3)?;
    }
    Ok(())
}

fn evaluate(
    fields: &Fields,
    p: SlenderJetParameters,
    limits: SlenderJetLimits,
    dx: f64,
    spatial: &mut Spatial,
    rhs: &mut Fields,
    gate: &CancelGate,
) -> Result<Evaluation, NonlinearJetError> {
    evaluate_split(fields, p, limits, dx, spatial, rhs, true, gate)
}

fn evaluate_split(
    fields: &Fields,
    p: SlenderJetParameters,
    limits: SlenderJetLimits,
    dx: f64,
    spatial: &mut Spatial,
    rhs: &mut Fields,
    explicit_viscosity: bool,
    gate: &CancelGate,
) -> Result<Evaluation, NonlinearJetError> {
    let n = fields.a.len();
    let mut volume = Sum::default();
    let mut momentum = Sum::default();
    let mut kinetic = Sum::default();
    let mut surface = Sum::default();
    let mut viscous = Sum::default();
    let mut upwind = Sum::default();
    let (mut min_h, mut max_h, mut max_slope, mut max_u) =
        (f64::INFINITY, 0.0_f64, 0.0_f64, 0.0_f64);
    for i in 0..n {
        poll_at(gate, i)?;
        let a = fields.a[i];
        if !a.is_finite() || a <= 0.0 || !fields.q[i].is_finite() {
            return Err(NonlinearJetError::InvalidState(
                "nonpositive area or nonfinite state",
            ));
        }
        let h = a.sqrt();
        spatial.h[i] = h;
        min_h = min_h.min(h);
        max_h = max_h.max(h);
        if h < limits.minimum_radius_m {
            return Err(domain(
                SlenderJetLimitKind::MinimumRadius,
                i,
                h,
                limits.minimum_radius_m,
            ));
        }
        let cells = h / dx;
        if cells < limits.minimum_radius_cells {
            return Err(domain(
                SlenderJetLimitKind::RadiusResolution,
                i,
                cells,
                limits.minimum_radius_cells,
            ));
        }
        volume.add(a);
        momentum.add(fields.q[i]);
    }
    for i in 0..n {
        poll_at(gate, i)?;
        let j = next(i, n);
        let u = finite(fields.q[i] / average(fields.a[i], fields.a[j]))?;
        spatial.u[i] = u;
        max_u = max_u.max(u.abs());
        let slope = finite((spatial.h[j] - spatial.h[i]) / dx)?;
        max_slope = max_slope.max(slope.abs());
        if slope.abs() > limits.maximum_abs_slope {
            return Err(domain(
                SlenderJetLimitKind::MaximumSlope,
                i,
                slope.abs(),
                limits.maximum_abs_slope,
            ));
        }
        kinetic.add((0.5 * fields.q[i]) * u);
        surface.add(average(spatial.h[i], spatial.h[j]) * (1.0 + slope * slope).sqrt());
    }
    let coefficient = p.surface_tension_n_m / p.density_kg_m3;
    let (mut capillary_jacobian, mut drainage_rate) = (0.0_f64, 0.0_f64);
    for i in 0..n {
        poll_at(gate, i)?;
        let l = previous(i, n);
        let r = next(i, n);
        let h = spatial.h[i];
        let a = fields.a[i];
        let mut capillary_stress = 0.0;
        if coefficient > 0.0 {
            let slope = (0.5 * (spatial.h[r] - spatial.h[l])) / dx;
            let second = ((spatial.h[r] - h) - (h - spatial.h[l])) / dx / dx;
            let inverse = 1.0 / (1.0 + slope * slope).sqrt();
            let inverse3 = inverse * inverse * inverse;
            let b = a * inverse3;
            // A*kappa - 2*h/sqrt(1+h_z^2), without canceling two h terms.
            capillary_stress = -h * inverse - b * second;
            // Absolute row sum of d(stress_i)/d(A_j), j=i-1,i,i+1.
            // Together with ||D_div D_grad||_infinity<=4/dx^2 this
            // bounds the frozen capillary operator used for the timestep.
            let slope_derivative =
                h * slope * inverse3 + 3.0 * b * second * slope * inverse * inverse;
            let center = (-inverse - 2.0 * h * second * inverse3 + 2.0 * b / dx / dx) / (2.0 * h);
            let left = (-0.5 * slope_derivative / dx - b / dx / dx) / (2.0 * spatial.h[l]);
            let right = (0.5 * slope_derivative / dx - b / dx / dx) / (2.0 * spatial.h[r]);
            capillary_jacobian =
                capillary_jacobian.max(finite(center.abs() + left.abs() + right.abs())?);
        }
        // This dual-cell mass flux exactly matches the averaged cell mass
        // update. A uniform velocity therefore remains uniform under pure
        // advection; upwinding dissipates kinetic energy instead of creating it.
        let mass_flux = average(fields.q[l], fields.q[i]);
        let velocity_jump = spatial.u[i] - spatial.u[l];
        let limited_slope = if mass_flux >= 0.0 {
            minmod(spatial.u[l] - spatial.u[previous(l, n)], velocity_jump)
        } else {
            minmod(velocity_jump, spatial.u[r] - spatial.u[i])
        };
        let advected_u = if mass_flux >= 0.0 {
            spatial.u[l] + 0.5 * limited_slope
        } else {
            spatial.u[i] - 0.5 * limited_slope
        };
        let strain = velocity_jump / dx;
        let viscous_stress = 3.0 * p.kinematic_viscosity_m2_s * a * strain;
        spatial.flux[i] = finite(
            mass_flux * advected_u + coefficient * capillary_stress
                - if explicit_viscosity {
                    viscous_stress
                } else {
                    0.0
                },
        )?;
        viscous.add(viscous_stress * strain);
        // minmod shares the jump's sign and has no larger magnitude, so this
        // kinetic loss is nonnegative. Reconstruction stays between the
        // upwind velocity and their centered average; it adds no material nu.
        upwind.add((0.5 * mass_flux.abs()) * velocity_jump * (velocity_jump - limited_slope));
        rhs.a[i] = finite(-(fields.q[i] - fields.q[l]) / dx)?;
        drainage_rate = drainage_rate.max(finite(rhs.a[i].abs() / a)?);
    }
    for i in 0..n {
        poll_at(gate, i)?;
        rhs.q[i] = finite(-(spatial.flux[next(i, n)] - spatial.flux[i]) / dx)?;
    }
    let pi = core::f64::consts::PI;
    let volume_m3 = positive(pi * dx * volume.total())?;
    let axial_momentum_kg_m_s = finite(pi * p.density_kg_m3 * dx * momentum.total())?;
    let diagnostics = SlenderJetDiagnostics {
        volume_m3,
        axial_momentum_kg_m_s,
        mean_velocity_m_s: finite(momentum.total() / volume.total())?,
        minimum_radius_m: min_h,
        maximum_radius_m: max_h,
        maximum_abs_slope: max_slope,
        kinetic_energy_j: finite(pi * p.density_kg_m3 * dx * kinetic.total())?,
        surface_energy_j: finite(2.0 * pi * p.surface_tension_n_m * dx * surface.total())?,
    };
    let capillary_rate = 2.0 * coefficient.sqrt() * capillary_jacobian.sqrt() / dx;
    let viscous_rate = if explicit_viscosity {
        12.0 * p.kinematic_viscosity_m2_s / dx / dx
    } else {
        0.0
    };
    let rate_per_s = finite(2.0 * max_u / dx + viscous_rate + capillary_rate + drainage_rate)?;
    Ok(Evaluation {
        diagnostics,
        rate_per_s,
        viscous_power_w: finite(pi * p.density_kg_m3 * dx * viscous.total())?,
        upwind_power_w: finite(pi * p.density_kg_m3 * upwind.total())?,
    })
}

fn blend(
    stage: &mut Fields,
    base: &Fields,
    rhs: &Fields,
    dt: f64,
    factor: f64,
    gate: &CancelGate,
) -> Result<(), NonlinearJetError> {
    for i in 0..stage.a.len() {
        poll_at(gate, i)?;
        // Difference form preserves an exactly stationary state.
        stage.a[i] = base.a[i] + factor * ((stage.a[i] - base.a[i]) + dt * rhs.a[i]);
        stage.q[i] = base.q[i] + factor * ((stage.q[i] - base.q[i]) + dt * rhs.q[i]);
    }
    Ok(())
}
fn stable_step(rate: f64, cfl: f64) -> f64 {
    if rate == 0.0 {
        f64::INFINITY
    } else {
        cfl / rate
    }
}
fn previous(i: usize, n: usize) -> usize {
    if i == 0 { n - 1 } else { i - 1 }
}
fn next(i: usize, n: usize) -> usize {
    if i + 1 == n { 0 } else { i + 1 }
}
fn minmod(a: f64, b: f64) -> f64 {
    if a > 0.0 && b > 0.0 {
        a.min(b)
    } else if a < 0.0 && b < 0.0 {
        a.max(b)
    } else {
        0.0
    }
}
fn average(a: f64, b: f64) -> f64 {
    if a == b { a } else { 0.5 * a + 0.5 * b }
}
fn finite(x: f64) -> Result<f64, NonlinearJetError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(NonlinearJetError::NonRepresentable)
    }
}
fn positive(x: f64) -> Result<f64, NonlinearJetError> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(NonlinearJetError::NonRepresentable)
    }
}
fn domain(kind: SlenderJetLimitKind, cell: usize, observed: f64, limit: f64) -> NonlinearJetError {
    NonlinearJetError::InitialDomain(SlenderJetViolation {
        kind,
        cell,
        observed,
        limit,
    })
}
fn poll(gate: &CancelGate) -> Result<(), NonlinearJetError> {
    if gate.is_requested() {
        Err(NonlinearJetError::Cancelled)
    } else {
        Ok(())
    }
}
fn poll_at(gate: &CancelGate, i: usize) -> Result<(), NonlinearJetError> {
    if i % 256 == 0 { poll(gate) } else { Ok(()) }
}
fn zeros(n: usize, gate: &CancelGate) -> Result<Vec<f64>, NonlinearJetError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| NonlinearJetError::Allocation)?;
    while values.len() < n {
        poll(gate)?;
        values.resize(n.min(values.len() + 256), 0.0);
    }
    Ok(values)
}
#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let corrected = value - self.correction;
        let next = self.sum + corrected;
        self.correction = (next - self.sum) - corrected;
        self.sum = next;
    }
    fn total(&self) -> f64 {
        self.sum
    }
}

#[cfg(test)]
mod tests;
