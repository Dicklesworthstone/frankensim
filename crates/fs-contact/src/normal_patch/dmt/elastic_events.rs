//! Sourced adhesive stick/slip driven by a finite tangential loading spring.
//!
//! Held original DMT area supplies the original constant-junction capacities.
//! Sliding obeys m*x''=Q-k*(x-r)-c*x'-F, with linear spring-stage position r.
//! The dashpot is tied to the fixed reference, so its loss is c*integral(v²).
//! Static support balances the actual spring load until its sourced threshold.
//! The loading spring is external apparatus, not invented junction compliance.
//!
//! Phase-limited entire-series cells evaluate the forced damped oscillator in
//! every damping regime. Each cell's acceleration has at most one zero; those
//! extrema partition velocity into monotone brackets, preventing skipped
//! underdamped stops. Only a bracketed event can set velocity to zero. Root,
//! phase, event, source and work limits refuse without mutating the checkpoint.
//! Numerical uncertainty is reported separately from genuine dissipative work;
//! this is not an interval certificate for nearly multiple event roots.

use super::events::{
    DmtSlipEventBudget, DmtSlipEventKind, DmtSlipEvents, DmtSlipRegime, MAX_SLIP_EVENTS,
};
use super::{DmtError, DmtResponse, SourceDmtSphereContact};
use crate::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_math::det;

#[path = "elastic_events/flow.rs"]
mod flow;
use flow::Flow;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtElasticEventBudget {
    pub motion: DmtSlipEventBudget,
    pub maximum_spring_extension_m: f64,
    pub maximum_stage_speed_m_s: f64,
    /// Includes phase-resolution cells and static segments; at most 4096.
    pub maximum_segments: usize,
    /// Per monotone root bracket; at most 96.
    pub maximum_root_iterations: usize,
    pub event_time_absolute_tolerance_s: f64,
    pub maximum_numerical_work_uncertainty_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtElasticEventState {
    identity: ContentHash,
    indentation_m: f64,
    tangent_displacement_m: f64,
    tangent_velocity_m_s: f64,
    stage_position_m: f64,
    time_s: f64,
    steps: u64,
}
impl DmtElasticEventState {
    pub const fn indentation_m(self) -> f64 {
        self.indentation_m
    }
    pub const fn tangent_displacement_m(self) -> f64 {
        self.tangent_displacement_m
    }
    pub const fn tangent_velocity_m_s(self) -> f64 {
        self.tangent_velocity_m_s
    }
    pub const fn stage_position_m(self) -> f64 {
        self.stage_position_m
    }
    pub const fn time_s(self) -> f64 {
        self.time_s
    }
    pub const fn steps(self) -> u64 {
        self.steps
    }
    pub fn regime(self) -> DmtSlipRegime {
        if self.indentation_m <= 0.0 {
            DmtSlipRegime::Released
        } else if self.tangent_velocity_m_s == 0.0 {
            DmtSlipRegime::Sticking
        } else {
            DmtSlipRegime::Sliding
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum DmtElasticEventError {
    Source(super::events::DmtSlipEventError),
    Invalid(&'static str),
    StateMismatch,
    Cancelled,
    Budget(&'static str),
    UnresolvedEvent,
    Balance {
        what: &'static str,
        residual: f64,
        tolerance: f64,
    },
}
impl core::fmt::Display for DmtElasticEventError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sourced elastic adhesive event refused: {self:?}")
    }
}
impl std::error::Error for DmtElasticEventError {}
impl From<super::events::DmtSlipEventError> for DmtElasticEventError {
    fn from(error: super::events::DmtSlipEventError) -> Self {
        Self::Source(error)
    }
}
impl From<DmtError> for DmtElasticEventError {
    fn from(error: DmtError) -> Self {
        Self::Source(error.into())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtElasticEvent {
    pub offset_s: f64,
    pub tangent_increment_m: f64,
    pub kind: DmtSlipEventKind,
    pub time_uncertainty_s: f64,
    pub stop_velocity_residual_m_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtElasticEventStep {
    state: DmtElasticEventState,
    pub contact: DmtResponse,
    pub tangent_increment_m: f64,
    pub accumulated_slip_m: f64,
    pub friction_impulse_n_s: f64,
    pub friction_time_moment_n_s2: f64,
    /// Signed resisting spring + fixed-reference damper + junction impulse.
    pub reaction_impulse_n_s: f64,
    pub reaction_time_moment_n_s2: f64,
    pub friction_dissipation_j: f64,
    pub damping_dissipation_j: f64,
    pub applied_tangent_work_j: f64,
    /// Work supplied at the moving spring anchor, -integral(k*(x-r)*r').
    pub stage_work_j: f64,
    /// Normal indentation is held throughout this lane, hence exactly zero.
    pub prescribed_normal_work_j: f64,
    pub spring_energy_j: f64,
    pub tangent_kinetic_energy_j: f64,
    pub total_energy_j: f64,
    pub energy_residual_j: f64,
    pub impulse_residual_n_s: f64,
    pub position_residual_m: f64,
    pub numerical_work_uncertainty_j: f64,
    pub events: [Option<DmtElasticEvent>; MAX_SLIP_EVENTS],
    pub event_count: usize,
    pub segments: usize,
    pub root_iterations: usize,
}
impl DmtElasticEventStep {
    pub const fn state(&self) -> DmtElasticEventState {
        self.state
    }
}

pub struct DmtElasticSlipEvents<'a> {
    source: DmtSlipEvents<'a>,
    stiffness_n_m: f64,
    damping_n_s_m: f64,
    budget: DmtElasticEventBudget,
    identity: ContentHash,
}
impl<'a> DmtElasticSlipEvents<'a> {
    pub fn new(
        normal: &'a SourceDmtSphereContact,
        shear: &'a SourceAdhesiveShear,
        mass_kg: f64,
        stiffness_n_m: f64,
        damping_n_s_m: f64,
        budget: DmtElasticEventBudget,
    ) -> Result<Self, DmtElasticEventError> {
        let source = DmtSlipEvents::new(normal, shear, mass_kg, budget.motion)?;
        for x in [
            stiffness_n_m,
            budget.maximum_spring_extension_m,
            budget.maximum_stage_speed_m_s,
            budget.event_time_absolute_tolerance_s,
            budget.maximum_numerical_work_uncertainty_j,
        ] {
            if !x.is_finite() || x <= 0.0 {
                return Err(DmtElasticEventError::Invalid(
                    "positive finite spring/budget",
                ));
            }
        }
        if !damping_n_s_m.is_finite()
            || damping_n_s_m < 0.0
            || budget.maximum_segments == 0
            || budget.maximum_segments > 4096
            || budget.maximum_root_iterations == 0
            || budget.maximum_root_iterations > 96
        {
            return Err(DmtElasticEventError::Invalid("damping/phase/root budget"));
        }
        finite(stiffness_n_m / mass_kg)?;
        finite(damping_n_s_m / mass_kg)?;
        let mut hash = DomainHasher::new("fs-contact/dmt-elastic-slip-events/v1");
        hash.update(normal.identity().as_bytes());
        hash.update(shear.resolved().identity().as_bytes());
        for x in [mass_kg, stiffness_n_m, damping_n_s_m] {
            hash.update(&x.to_bits().to_le_bytes());
        }
        Ok(Self {
            source,
            stiffness_n_m,
            damping_n_s_m,
            budget,
            identity: hash.finalize(),
        })
    }
    pub const fn normal(&self) -> &'a SourceDmtSphereContact {
        self.source.normal()
    }
    pub const fn shear(&self) -> &'a SourceAdhesiveShear {
        self.source.shear()
    }
    pub const fn mass_kg(&self) -> f64 {
        self.source.mass_kg()
    }
    pub const fn stiffness_n_per_m(&self) -> f64 {
        self.stiffness_n_m
    }
    pub const fn damping_n_s_per_m(&self) -> f64 {
        self.damping_n_s_m
    }
    pub const fn budget(&self) -> DmtElasticEventBudget {
        self.budget
    }

    pub fn initialize(
        &self,
        indentation_m: f64,
        tangent_displacement_m: f64,
        tangent_velocity_m_s: f64,
        stage_position_m: f64,
    ) -> Result<DmtElasticEventState, DmtElasticEventError> {
        self.source
            .initialize(indentation_m, tangent_displacement_m, tangent_velocity_m_s)?;
        self.check_extension(finite(tangent_displacement_m - stage_position_m)?)?;
        self.check_force(self.damping_n_s_m * tangent_velocity_m_s)?;
        Ok(DmtElasticEventState {
            identity: self.identity,
            indentation_m,
            tangent_displacement_m,
            tangent_velocity_m_s,
            stage_position_m,
            time_s: 0.0,
            steps: 0,
        })
    }

    /// Held indentation and applied force; the spring anchor moves linearly.
    /// Retain only the returned checkpoint on success. No runtime allocation.
    pub fn advance(
        &self,
        previous: DmtElasticEventState,
        next_stage_m: f64,
        applied_tangent_force_n: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<DmtElasticEventStep, DmtElasticEventError> {
        poll(gate)?;
        if previous.identity != self.identity {
            return Err(DmtElasticEventError::StateMismatch);
        }
        if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > self.budget.motion.maximum_step_s {
            return Err(DmtElasticEventError::Budget("step duration"));
        }
        self.source.initialize(
            previous.indentation_m,
            previous.tangent_displacement_m,
            previous.tangent_velocity_m_s,
        )?;
        self.check_force(applied_tangent_force_n)?;
        let u = finite((next_stage_m - previous.stage_position_m) / dt_s)?;
        bounded(
            u.abs(),
            self.budget.maximum_stage_speed_m_s,
            "spring stage speed",
        )?;
        let contact = self.normal().evaluate(previous.indentation_m)?;
        let capacity = self.shear().law().static_shear_strength_pa() * contact.contact_area_m2;
        let kinetic = self.shear().law().kinetic_shear_strength_pa() * contact.contact_area_m2;
        let mut result = DmtElasticEventStep {
            state: previous,
            contact,
            tangent_increment_m: 0.0,
            accumulated_slip_m: 0.0,
            friction_impulse_n_s: 0.0,
            friction_time_moment_n_s2: 0.0,
            reaction_impulse_n_s: 0.0,
            reaction_time_moment_n_s2: 0.0,
            friction_dissipation_j: 0.0,
            damping_dissipation_j: 0.0,
            applied_tangent_work_j: 0.0,
            stage_work_j: 0.0,
            prescribed_normal_work_j: 0.0,
            spring_energy_j: 0.0,
            tangent_kinetic_energy_j: 0.0,
            total_energy_j: 0.0,
            energy_residual_j: 0.0,
            impulse_residual_n_s: 0.0,
            position_residual_m: 0.0,
            numerical_work_uncertainty_j: 0.0,
            events: [None; MAX_SLIP_EVENTS],
            event_count: 0,
            segments: 0,
            root_iterations: 0,
        };
        let m = self.mass_kg();
        let k = self.stiffness_n_m;
        let c = self.damping_n_s_m;
        let omega2 = finite(k / m)?;
        let drag = finite(c / m)?;
        let phase_step = finite(0.25 / finite(det::sqrt(omega2) + drag)?)?;
        let initial_extension =
            finite(previous.tangent_displacement_m - previous.stage_position_m)?;
        self.check_extension(initial_extension)?;
        let initial_energy = finite(
            0.5 * m * previous.tangent_velocity_m_s.powi(2) + 0.5 * k * initial_extension.powi(2),
        )?;
        let mut v = previous.tangent_velocity_m_s;
        let mut offset = 0.0;
        while offset < dt_s {
            poll(gate)?;
            if result.segments >= self.budget.maximum_segments {
                return Err(DmtElasticEventError::Budget("phase segments"));
            }
            let x = finite(previous.tangent_displacement_m + result.tangent_increment_m)?;
            let r = finite(previous.stage_position_m + u * offset)?;
            let extension = finite(x - r)?;
            self.check_extension(extension)?;
            let load = finite(applied_tangent_force_n - k * extension)?;
            let remaining = dt_s - offset;
            let mut direction = v.signum();
            if v == 0.0 {
                let slope = finite(k * u)?;
                let outward = load.abs() == capacity && load * slope > 0.0;
                let zero_threshold = capacity == 0.0 && load == 0.0 && slope != 0.0;
                if load.abs() <= capacity && !outward && !zero_threshold {
                    let release = if slope > 0.0 {
                        (capacity - load) / slope
                    } else if slope < 0.0 {
                        (-capacity - load) / slope
                    } else {
                        f64::INFINITY
                    };
                    let h = remaining.min(release);
                    if h <= 0.0 || offset + h <= offset {
                        return Err(DmtElasticEventError::UnresolvedEvent);
                    }
                    self.check_extension(extension - u * h)?;
                    let js = finite(k * h * (extension - 0.5 * u * h))?;
                    let jf = finite(h * (load + 0.5 * slope * h))?;
                    let mf = finite(h * h * (0.5 * load + slope * h / 3.0))?;
                    let end_extension = extension - u * h;
                    let coordinate_error =
                        64.0 * f64::EPSILON * (x.abs() + r.abs() + (u * h).abs());
                    let static_uncertainty = finite(
                        128.0
                            * f64::EPSILON
                            * ((u * js).abs()
                                + 0.5 * k * (extension.powi(2) + end_extension.powi(2)))
                            + (k * (extension.abs() + end_extension.abs())
                                + applied_tangent_force_n.abs())
                                * coordinate_error,
                    )?;
                    self.accumulate(
                        &mut result,
                        offset,
                        0.0,
                        0.0,
                        jf,
                        mf,
                        applied_tangent_force_n * h,
                        0.5 * applied_tangent_force_n * h * h,
                        0.0,
                        0.0,
                        -u * js,
                        static_uncertainty,
                    )?;
                    offset += h;
                    if h == remaining {
                        break;
                    }
                    direction = slope.signum();
                } else {
                    direction = if load == 0.0 {
                        slope.signum()
                    } else {
                        load.signum()
                    };
                }
                self.event(&mut result, offset, DmtSlipEventKind::SlipStarted, 0.0, 0.0)?;
            }
            if offset >= dt_s {
                break;
            }
            let x = finite(previous.tangent_displacement_m + result.tangent_increment_m)?;
            let extension = finite(x - (previous.stage_position_m + u * offset))?;
            let friction = direction * kinetic;
            let acceleration =
                finite((applied_tangent_force_n - k * extension - c * v - friction) / m)?;
            let h_max = (dt_s - offset).min(phase_step);
            if h_max <= 0.0 || offset + h_max <= offset {
                return Err(DmtElasticEventError::UnresolvedEvent);
            }
            let flow = Flow::new(v, acceleration, u, omega2, drag, h_max)?;
            let extremum = self.extremum(|z| flow.a(z), 1.0, h_max, &mut result, gate)?;
            let stop = self.stop(&flow, direction, extremum, &mut result, gate)?;
            let (z, time_error) = stop.unwrap_or((1.0, 0.0));
            let h = h_max * z;
            if h <= 0.0 || offset + h <= offset {
                return Err(DmtElasticEventError::UnresolvedEvent);
            }
            self.check_path(
                &flow,
                z,
                extremum,
                extension,
                u,
                friction,
                applied_tangent_force_n,
                &mut result,
                gate,
            )?;
            let dx = finite(flow.dx(z))?;
            let slip = finite(direction * dx)?;
            if slip < 0.0 {
                return Err(DmtElasticEventError::UnresolvedEvent);
            }
            let [tv, integral_x, integral_tx] = flow.moments(z);
            let js = finite(k * (extension * h + integral_x - 0.5 * u * h * h))?;
            let ms = finite(k * (0.5 * extension * h * h + integral_tx - u * h * h * h / 3.0))?;
            let jf = finite(friction * h)?;
            let mf = finite(0.5 * friction * h * h)?;
            let jr = finite(js + c * dx + jf)?;
            let mr = finite(ms + c * tv + mf)?;
            let (damping, uncertainty) = flow.damping_loss(c, z, u)?;
            let loss = finite(kinetic * slip)?;
            let stage_work = finite(-u * js)?;
            let endpoint = finite(flow.v(z))?;
            let end_extension = extension + dx - u * h;
            // Propagation uncertainty is present even in a conservative cell.
            // Estimate series/force arithmetic error, then carry displacement,
            // impulse and time-moment errors to kinetic/spring/work scales.
            // This is a numerical estimate, not an interval certificate.
            let round = 64.0 * f64::EPSILON;
            let force_scale =
                applied_tangent_force_n.abs() + (k * extension).abs() + (c * v).abs() + kinetic;
            let velocity_error = 3.0 * (flow.speed_roundoff() + round * h_max / m * force_scale);
            let dx_error = h * velocity_error
                + round
                    * (x.abs()
                        + (previous.stage_position_m + u * offset).abs()
                        + dx.abs()
                        + (u * h).abs());
            let impulse_error = k * h * dx_error
                + c * h * velocity_error
                + round * (js.abs() + (c * dx).abs() + jf.abs());
            let moment_error = h * impulse_error + round * (ms.abs() + (c * tv).abs() + mf.abs());
            let position_error = dx_error + (h * impulse_error + moment_error) / m;
            let speed_error = velocity_error + impulse_error / m;
            let speed_scale = flow.speed_bound() + flow.speed_roundoff();
            let energy_scale = 0.5 * m * (v * v + endpoint * endpoint)
                + 0.5 * k * (extension * extension + end_extension * end_extension)
                + stage_work.abs()
                + loss
                + damping
                + applied_tangent_force_n.abs() * slip;
            let uncertainty = finite(
                uncertainty
                    + m * (v.abs() + endpoint.abs() + speed_error) * speed_error
                    + k * (extension.abs() + end_extension.abs() + position_error) * position_error
                    + (applied_tangent_force_n.abs() + kinetic) * position_error
                    + u.abs() * k * h * dx_error
                    + c * h * (2.0 * speed_scale + speed_error) * speed_error
                    + 2.0 * round * energy_scale,
            )?;
            self.accumulate(
                &mut result,
                offset,
                dx,
                slip,
                jf,
                mf,
                jr,
                mr,
                loss,
                damping,
                stage_work,
                uncertainty,
            )?;
            offset += h;
            if stop.is_some() {
                let tolerance = self.impulse_tolerance(
                    m * v.abs() + applied_tangent_force_n.abs() * h + jr.abs(),
                )?;
                balance("event velocity", m * endpoint, tolerance)?;
                result.numerical_work_uncertainty_j =
                    finite(result.numerical_work_uncertainty_j + 0.5 * m * endpoint * endpoint)?;
                v = 0.0;
                self.event(
                    &mut result,
                    offset,
                    DmtSlipEventKind::Stopped,
                    time_error,
                    endpoint,
                )?;
            } else {
                v = endpoint;
            }
        }
        let end_x = finite(previous.tangent_displacement_m + result.tangent_increment_m)?;
        let extension = finite(end_x - next_stage_m)?;
        self.check_extension(extension)?;
        result.spring_energy_j = finite(0.5 * k * extension * extension)?;
        result.tangent_kinetic_energy_j = finite(0.5 * m * v * v)?;
        let after = finite(result.spring_energy_j + result.tangent_kinetic_energy_j)?;
        result.total_energy_j = finite(after + contact.reversible_energy_j)?;
        result.applied_tangent_work_j =
            finite(applied_tangent_force_n * result.tangent_increment_m)?;
        result.energy_residual_j = finite(
            after - initial_energy + result.friction_dissipation_j + result.damping_dissipation_j
                - result.applied_tangent_work_j
                - result.stage_work_j,
        )?;
        let energy_scale = initial_energy.abs()
            + after.abs()
            + result.friction_dissipation_j
            + result.damping_dissipation_j
            + result.applied_tangent_work_j.abs()
            + result.stage_work_j.abs();
        let energy_tolerance = finite(
            self.budget.motion.energy_absolute_tolerance_j
                + self.budget.motion.energy_relative_tolerance * energy_scale,
        )?;
        bounded(
            result.numerical_work_uncertainty_j,
            self.budget.maximum_numerical_work_uncertainty_j,
            "numerical work uncertainty",
        )?;
        balance(
            "energy plus numerical uncertainty",
            result.energy_residual_j.abs() + result.numerical_work_uncertainty_j,
            energy_tolerance,
        )?;
        let momentum = finite(m * (v - previous.tangent_velocity_m_s))?;
        result.impulse_residual_n_s =
            finite(momentum - applied_tangent_force_n * dt_s + result.reaction_impulse_n_s)?;
        let impulse_tolerance = self.impulse_tolerance(
            momentum.abs()
                + (applied_tangent_force_n * dt_s).abs()
                + result.reaction_impulse_n_s.abs(),
        )?;
        balance(
            "reaction impulse",
            result.impulse_residual_n_s,
            impulse_tolerance,
        )?;
        result.position_residual_m = finite(
            result.tangent_increment_m
                - dt_s * previous.tangent_velocity_m_s
                - 0.5 * applied_tangent_force_n / m * dt_s * dt_s
                + (dt_s * result.reaction_impulse_n_s - result.reaction_time_moment_n_s2) / m,
        )?;
        balance(
            "reaction first moment",
            result.position_residual_m,
            impulse_tolerance * dt_s / m,
        )?;
        result.state = DmtElasticEventState {
            identity: self.identity,
            indentation_m: previous.indentation_m,
            tangent_displacement_m: end_x,
            tangent_velocity_m_s: v,
            stage_position_m: next_stage_m,
            time_s: finite(previous.time_s + dt_s)?,
            steps: previous
                .steps
                .checked_add(1)
                .ok_or(DmtElasticEventError::Budget("step counter"))?,
        };
        if result.state.time_s <= previous.time_s {
            return Err(DmtElasticEventError::Budget("clock resolution"));
        }
        poll(gate)?;
        Ok(result)
    }

    fn impulse_tolerance(&self, scale: f64) -> Result<f64, DmtElasticEventError> {
        finite(
            self.budget.motion.impulse_absolute_tolerance_n_s
                + self.budget.motion.impulse_relative_tolerance * scale,
        )
    }
    fn root(
        &self,
        value: impl Fn(f64) -> f64,
        mut lo: f64,
        mut hi: f64,
        h: f64,
        result: &mut DmtElasticEventStep,
        gate: &CancelGate,
    ) -> Result<(f64, f64), DmtElasticEventError> {
        let mut a = finite(value(lo))?;
        let b = finite(value(hi))?;
        if a == 0.0 {
            return Ok((lo, 0.0));
        }
        if b == 0.0 {
            return Ok((hi, 0.0));
        }
        if a.signum() == b.signum() {
            return Err(DmtElasticEventError::UnresolvedEvent);
        }
        for _ in 0..self.budget.maximum_root_iterations {
            poll(gate)?;
            result.root_iterations += 1;
            let mid = f64::midpoint(lo, hi);
            let v = finite(value(mid))?;
            if v == 0.0 {
                return Ok((mid, 0.0));
            }
            if (hi - lo) * h <= self.budget.event_time_absolute_tolerance_s {
                return Ok((mid, (hi - lo) * h));
            }
            if mid == lo || mid == hi {
                return Err(DmtElasticEventError::UnresolvedEvent);
            }
            if v.signum() == a.signum() {
                lo = mid;
                a = v;
            } else {
                hi = mid;
            }
        }
        Err(DmtElasticEventError::Budget("event root iterations"))
    }
    fn extremum(
        &self,
        derivative: impl Fn(f64) -> f64,
        end: f64,
        h: f64,
        result: &mut DmtElasticEventStep,
        gate: &CancelGate,
    ) -> Result<Option<f64>, DmtElasticEventError> {
        let a = finite(derivative(0.0))?;
        let b = finite(derivative(end))?;
        if a != 0.0 && b != 0.0 && a.signum() != b.signum() {
            Ok(Some(self.root(derivative, 0.0, end, h, result, gate)?.0))
        } else {
            Ok(None)
        }
    }
    fn stop(
        &self,
        flow: &Flow,
        direction: f64,
        extremum: Option<f64>,
        result: &mut DmtElasticEventStep,
        gate: &CancelGate,
    ) -> Result<Option<(f64, f64)>, DmtElasticEventError> {
        let mut lo = 0.0;
        for hi in [extremum.unwrap_or(1.0), 1.0] {
            if hi <= lo {
                continue;
            }
            let a = finite(direction * flow.v(lo))?;
            let b = finite(direction * flow.v(hi))?;
            if hi < 1.0 && b != 0.0 && b.abs() <= flow.speed_roundoff() {
                return Err(DmtElasticEventError::UnresolvedEvent);
            }
            if a > 0.0 && b <= 0.0 {
                return Ok(Some(self.root(
                    |z| flow.v(z),
                    lo,
                    hi,
                    flow.h,
                    result,
                    gate,
                )?));
            }
            if a < 0.0 || b < 0.0 {
                return Err(DmtElasticEventError::UnresolvedEvent);
            }
            if a == 0.0 && lo > 0.0 {
                return Ok(Some((lo, self.budget.event_time_absolute_tolerance_s)));
            }
            lo = hi;
        }
        Ok(None)
    }
    #[allow(clippy::too_many_arguments)]
    fn check_path(
        &self,
        flow: &Flow,
        end: f64,
        extremum: Option<f64>,
        extension: f64,
        u: f64,
        friction: f64,
        applied: f64,
        result: &mut DmtElasticEventStep,
        gate: &CancelGate,
    ) -> Result<(), DmtElasticEventError> {
        let check = |z: f64| -> Result<(), DmtElasticEventError> {
            let v = finite(flow.v(z))?;
            bounded(
                v.abs(),
                self.budget.motion.maximum_tangent_speed_m_s,
                "tangent speed",
            )?;
            self.check_extension(extension + flow.dx(z) - u * flow.h * z)?;
            self.check_force(self.damping_n_s_m * v)?;
            self.check_force(applied - self.mass_kg() * flow.a(z) / flow.h)?;
            self.check_force(friction)
        };
        check(0.0)?;
        check(end)?;
        let peak = extremum.filter(|z| *z < end);
        if let Some(z) = peak {
            check(z)?;
        }
        if let Some(z) = self.extremum(|z| flow.jerk(z), end, flow.h, result, gate)? {
            check(z)?;
        }
        // Spring extension extrema satisfy v=u. Its roots are found on the
        // same monotone velocity pieces, not guessed from endpoint loads.
        let mut lo = 0.0;
        for hi in [peak.unwrap_or(end), end] {
            if hi <= lo {
                continue;
            }
            let a = flow.v(lo) - u;
            let b = flow.v(hi) - u;
            if a != 0.0 && b != 0.0 && a.signum() != b.signum() {
                check(
                    self.root(|z| flow.v(z) - u, lo, hi, flow.h, result, gate)?
                        .0,
                )?;
            }
            lo = hi;
        }
        Ok(())
    }
    #[allow(clippy::too_many_arguments)]
    fn accumulate(
        &self,
        r: &mut DmtElasticEventStep,
        offset: f64,
        dx: f64,
        slip: f64,
        jf: f64,
        mf: f64,
        jr: f64,
        mr: f64,
        friction: f64,
        damping: f64,
        stage: f64,
        uncertainty: f64,
    ) -> Result<(), DmtElasticEventError> {
        if r.segments >= self.budget.maximum_segments {
            return Err(DmtElasticEventError::Budget("phase segments"));
        }
        r.tangent_increment_m = finite(r.tangent_increment_m + dx)?;
        r.accumulated_slip_m = finite(r.accumulated_slip_m + slip)?;
        r.friction_impulse_n_s = finite(r.friction_impulse_n_s + jf)?;
        r.friction_time_moment_n_s2 = finite(r.friction_time_moment_n_s2 + offset * jf + mf)?;
        r.reaction_impulse_n_s = finite(r.reaction_impulse_n_s + jr)?;
        r.reaction_time_moment_n_s2 = finite(r.reaction_time_moment_n_s2 + offset * jr + mr)?;
        r.friction_dissipation_j = finite(r.friction_dissipation_j + friction)?;
        r.damping_dissipation_j = finite(r.damping_dissipation_j + damping)?;
        r.stage_work_j = finite(r.stage_work_j + stage)?;
        r.numerical_work_uncertainty_j = finite(r.numerical_work_uncertainty_j + uncertainty)?;
        r.segments += 1;
        bounded(
            r.accumulated_slip_m,
            self.budget.motion.maximum_accumulated_slip_m,
            "total absolute slip",
        )
    }
    fn event(
        &self,
        r: &mut DmtElasticEventStep,
        time: f64,
        kind: DmtSlipEventKind,
        uncertainty: f64,
        residual: f64,
    ) -> Result<(), DmtElasticEventError> {
        if r.event_count >= self.budget.motion.maximum_events {
            return Err(DmtElasticEventError::Budget("events"));
        }
        r.events[r.event_count] = Some(DmtElasticEvent {
            offset_s: time,
            tangent_increment_m: r.tangent_increment_m,
            kind,
            time_uncertainty_s: uncertainty,
            stop_velocity_residual_m_s: residual,
        });
        r.event_count += 1;
        Ok(())
    }
    fn check_extension(&self, x: f64) -> Result<(), DmtElasticEventError> {
        bounded(
            x.abs(),
            self.budget.maximum_spring_extension_m,
            "spring extension",
        )?;
        self.check_force(self.stiffness_n_m * x)
    }
    fn check_force(&self, force: f64) -> Result<(), DmtElasticEventError> {
        bounded(force.abs(), self.budget.motion.maximum_force_n, "force")
    }
}

fn poll(gate: &CancelGate) -> Result<(), DmtElasticEventError> {
    if gate.is_requested() {
        Err(DmtElasticEventError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, DmtElasticEventError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DmtElasticEventError::Invalid("finite elastic trajectory"))
    }
}
fn bounded(value: f64, ceiling: f64, what: &'static str) -> Result<(), DmtElasticEventError> {
    finite(value)?;
    if value <= ceiling {
        Ok(())
    } else {
        Err(DmtElasticEventError::Budget(what))
    }
}
fn balance(what: &'static str, residual: f64, tolerance: f64) -> Result<(), DmtElasticEventError> {
    finite(residual)?;
    finite(tolerance)?;
    if residual.abs() <= tolerance {
        Ok(())
    } else {
        Err(DmtElasticEventError::Balance {
            what,
            residual,
            tolerance,
        })
    }
}
