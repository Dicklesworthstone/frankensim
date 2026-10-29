//! Event-resolved collinear junction sliding with prescribed DMT normal motion.
//!
//! The two surfaces have a fixed tangent direction and finite relative mass.
//! Signed normal approach is linear in time over a call, while the applied
//! tangential force is held. Original DMT geometry gives A=pi R max(delta,0).
//! Consequently the original constant-junction strengths give piecewise affine
//! static/kinetic capacities. Sliding velocity is quadratic in time: actual
//! stops, static release, and geometric contact entry/exit can be located
//! without a velocity capture threshold or a timestep-dependent force clamp.
//!
//! Each sign-preserving sliding segment integrates impulse, its first time
//! moment, displacement and friction work analytically. Positive two-point
//! Gaussian integration is exact for the cubic work integrand. A stopped
//! segment supports the applied force only inside the current static cone.
//! No shear spring, pressure-friction term, heat recipient, JKR peeling,
//! arbitrary vector stopping rule or normal inertial trajectory is inferred.
//! Normal motion is an explicit ideal constraint; its supplied work equals
//! the original DMT reversible work. Floating-point gates are not interval
//! certificates for nearly multiple event roots.

use super::{DmtError, DmtResponse, SourceDmtSphereContact};
use crate::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_math::det;

pub const MAX_SLIP_EVENTS: usize = 16;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtSlipEventBudget {
    pub maximum_step_s: f64,
    pub minimum_indentation_m: f64,
    pub maximum_normal_increment_m: f64,
    /// Actual absolute path length, including motion before and after reversal.
    pub maximum_accumulated_slip_m: f64,
    pub maximum_tangent_speed_m_s: f64,
    /// Includes applied force, normal response, and both shear capacities.
    pub maximum_force_n: f64,
    pub maximum_events: usize,
    pub impulse_absolute_tolerance_n_s: f64,
    pub impulse_relative_tolerance: f64,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DmtSlipEventKind {
    Stopped,
    SlipStarted,
    ContactEntered,
    ContactLeft,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtSlipEvent {
    /// Time from the beginning of this call, s.
    pub offset_s: f64,
    pub tangent_increment_m: f64,
    pub kind: DmtSlipEventKind,
    /// Polynomial residual before publishing an analytically located stop.
    pub stop_velocity_residual_m_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DmtSlipRegime {
    Released,
    Sticking,
    Sliding,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtSlipEventState {
    identity: ContentHash,
    indentation_m: f64,
    tangent_displacement_m: f64,
    tangent_velocity_m_s: f64,
    time_s: f64,
    steps: u64,
}
impl DmtSlipEventState {
    pub const fn indentation_m(self) -> f64 {
        self.indentation_m
    }
    pub const fn tangent_displacement_m(self) -> f64 {
        self.tangent_displacement_m
    }
    pub const fn tangent_velocity_m_s(self) -> f64 {
        self.tangent_velocity_m_s
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
pub enum DmtSlipEventError {
    Normal(DmtError),
    Invalid(&'static str),
    StateMismatch,
    Cancelled,
    EventBudget,
    Range(&'static str),
    UnresolvedEvent,
    ImpulseBalance {
        residual_n_s: f64,
        tolerance_n_s: f64,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl From<DmtError> for DmtSlipEventError {
    fn from(value: DmtError) -> Self {
        Self::Normal(value)
    }
}
impl core::fmt::Display for DmtSlipEventError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sourced DMT slip event refused: {self:?}")
    }
}
impl std::error::Error for DmtSlipEventError {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtSlipEventStep {
    state: DmtSlipEventState,
    pub contact: DmtResponse,
    pub initial_contact_area_m2: f64,
    pub tangent_increment_m: f64,
    pub accumulated_slip_m: f64,
    /// Signed reaction impulse; the force on the relative mass is its negative.
    pub friction_impulse_n_s: f64,
    /// Integral of t F_friction(t) dt, with t measured from this call's start.
    pub friction_time_moment_n_s2: f64,
    pub friction_dissipation_j: f64,
    pub applied_tangent_work_j: f64,
    /// Supplied by the ideal normal-motion constraint, retaining signed work.
    pub prescribed_normal_work_j: f64,
    pub tangent_kinetic_energy_j: f64,
    pub total_energy_j: f64,
    pub energy_residual_j: f64,
    pub impulse_residual_n_s: f64,
    pub events: [Option<DmtSlipEvent>; MAX_SLIP_EVENTS],
    pub event_count: usize,
    pub segments: usize,
}
impl DmtSlipEventStep {
    pub const fn state(&self) -> DmtSlipEventState {
        self.state
    }
}

/// Fixed-source, fixed-tangent finite relative mass. A pair of unconstrained
/// bodies uses reduced mass m1*m2/(m1+m2); common-mode motion stays with the
/// receiving mechanics. Construction and stepping allocate nothing.
pub struct DmtSlipEvents<'a> {
    normal: &'a SourceDmtSphereContact,
    shear: &'a SourceAdhesiveShear,
    mass_kg: f64,
    budget: DmtSlipEventBudget,
    identity: ContentHash,
}

impl<'a> DmtSlipEvents<'a> {
    pub fn new(
        normal: &'a SourceDmtSphereContact,
        shear: &'a SourceAdhesiveShear,
        mass_kg: f64,
        budget: DmtSlipEventBudget,
    ) -> Result<Self, DmtSlipEventError> {
        if shear.normal_state_identity() != normal.normal_state().identity() {
            return Err(DmtSlipEventError::StateMismatch);
        }
        for v in [
            mass_kg,
            budget.maximum_step_s,
            budget.maximum_normal_increment_m,
            budget.maximum_accumulated_slip_m,
            budget.maximum_tangent_speed_m_s,
            budget.maximum_force_n,
            budget.impulse_absolute_tolerance_n_s,
            budget.energy_absolute_tolerance_j,
        ] {
            if !v.is_finite() || v <= 0.0 {
                return Err(DmtSlipEventError::Invalid("positive finite budget/mass"));
            }
        }
        for v in [
            budget.impulse_relative_tolerance,
            budget.energy_relative_tolerance,
        ] {
            if !v.is_finite() || v < 0.0 {
                return Err(DmtSlipEventError::Invalid(
                    "nonnegative finite relative tolerance",
                ));
            }
        }
        if !budget.minimum_indentation_m.is_finite()
            || budget.minimum_indentation_m >= 0.0
            || budget.maximum_events == 0
            || budget.maximum_events > MAX_SLIP_EVENTS
        {
            return Err(DmtSlipEventError::Invalid("normal/event budget"));
        }
        let mut hash = DomainHasher::new("fs-contact/dmt-collinear-slip-events/v1");
        hash.update(normal.identity().as_bytes());
        hash.update(shear.resolved().identity().as_bytes());
        hash.update(&mass_kg.to_bits().to_le_bytes());
        Ok(Self {
            normal,
            shear,
            mass_kg,
            budget,
            identity: hash.finalize(),
        })
    }
    pub const fn normal(&self) -> &'a SourceDmtSphereContact {
        self.normal
    }
    pub const fn shear(&self) -> &'a SourceAdhesiveShear {
        self.shear
    }
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    pub const fn budget(&self) -> DmtSlipEventBudget {
        self.budget
    }

    pub fn initialize(
        &self,
        indentation_m: f64,
        tangent_displacement_m: f64,
        tangent_velocity_m_s: f64,
    ) -> Result<DmtSlipEventState, DmtSlipEventError> {
        self.check_normal(indentation_m)?;
        self.check_velocity(tangent_velocity_m_s)?;
        finite(tangent_displacement_m)?;
        Ok(DmtSlipEventState {
            identity: self.identity,
            indentation_m,
            tangent_displacement_m,
            tangent_velocity_m_s,
            time_s: 0.0,
            steps: 0,
        })
    }

    /// Integrate the actual collinear trajectory for held applied tangent force
    /// and linearly prescribed normal approach. Accept by retaining the returned
    /// state; failed or cancelled proposals never change their input checkpoint.
    pub fn advance(
        &self,
        previous: DmtSlipEventState,
        next_indentation_m: f64,
        applied_tangent_force_n: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<DmtSlipEventStep, DmtSlipEventError> {
        poll(gate)?;
        if previous.identity != self.identity {
            return Err(DmtSlipEventError::StateMismatch);
        }
        if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > self.budget.maximum_step_s {
            return Err(DmtSlipEventError::Range("step duration"));
        }
        self.check_force(applied_tangent_force_n)?;
        let start = self.check_normal(previous.indentation_m)?;
        let contact = self.check_normal(next_indentation_m)?;
        let delta = finite(next_indentation_m - previous.indentation_m)?;
        if delta.abs() > self.budget.maximum_normal_increment_m {
            return Err(DmtSlipEventError::Range("normal increment"));
        }
        let normal_rate = finite(delta / dt_s)?;
        let time = finite(previous.time_s + dt_s)?;
        if time <= previous.time_s {
            return Err(DmtSlipEventError::Range("clock resolution"));
        }
        let steps = previous
            .steps
            .checked_add(1)
            .ok_or(DmtSlipEventError::Range("step counter"))?;
        // The source force is monotone on each side of geometric touch; its
        // largest attractive magnitude is at touch. Check that interior extremum.
        let crosses = (previous.indentation_m < 0.0 && next_indentation_m >= 0.0)
            || (previous.indentation_m > 0.0 && next_indentation_m <= 0.0);
        let boundary = if crosses {
            self.check_normal(0.0)?;
            Some(dt_s * (-previous.indentation_m / delta))
        } else {
            None
        };
        let mut result = DmtSlipEventStep {
            state: previous,
            contact,
            initial_contact_area_m2: start.contact_area_m2,
            tangent_increment_m: 0.0,
            accumulated_slip_m: 0.0,
            friction_impulse_n_s: 0.0,
            friction_time_moment_n_s2: 0.0,
            friction_dissipation_j: 0.0,
            applied_tangent_work_j: 0.0,
            prescribed_normal_work_j: self
                .normal
                .reversible_work_between(previous.indentation_m, next_indentation_m)?,
            tangent_kinetic_energy_j: 0.0,
            total_energy_j: 0.0,
            energy_residual_j: 0.0,
            impulse_residual_n_s: 0.0,
            events: [None; MAX_SLIP_EVENTS],
            event_count: 0,
            segments: 0,
        };
        let mut offset = 0.0;
        let mut velocity = previous.tangent_velocity_m_s;
        let mut contact_boundary_done = false;
        while offset < dt_s {
            poll(gate)?;
            // Every zero-duration transition consumes an event. Positive
            // segments either finish or end at one of those bounded events.
            if result.segments > 2 * self.budget.maximum_events + 3 {
                return Err(DmtSlipEventError::EventBudget);
            }
            let phase_end = if !contact_boundary_done {
                boundary.unwrap_or(dt_s)
            } else {
                dt_s
            };
            if phase_end == offset {
                self.push_event(
                    &mut result,
                    offset,
                    if normal_rate > 0.0 {
                        DmtSlipEventKind::ContactEntered
                    } else {
                        DmtSlipEventKind::ContactLeft
                    },
                    0.0,
                )?;
                contact_boundary_done = true;
                continue;
            }
            let h_max = phase_end - offset;
            if h_max <= 0.0 {
                return Err(DmtSlipEventError::UnresolvedEvent);
            }
            // The physical DMT area is affine on each side of touch. Anchor
            // that line to ORIGINAL admitted areas, rather than rearranging
            // pi*a^2 as pi*R*delta: an ulp difference at an exact static
            // threshold must not trigger a finite jump to kinetic strength.
            let (phase_start, area_start, area_end) = match boundary {
                Some(touch) if contact_boundary_done => (touch, 0.0, contact.contact_area_m2),
                Some(_) => (0.0, start.contact_area_m2, 0.0),
                None => (0.0, start.contact_area_m2, contact.contact_area_m2),
            };
            let area_rate = (area_end - area_start) / (phase_end - phase_start);
            let area = area_start + area_rate * (offset - phase_start);
            let ks = self.shear.law().static_shear_strength_pa();
            let kk = self.shear.law().kinetic_shear_strength_pa();
            let static_force = ks * area;
            let static_rate = ks * area_rate;
            let kinetic_force = kk * area;
            let kinetic_rate = kk * area_rate;
            let force = applied_tangent_force_n;
            self.check_force(static_force)?;
            self.check_force(kinetic_force)?;
            // Equality with a declining static cone releases immediately;
            // equality with a fixed or increasing cone remains supported.
            let held = velocity == 0.0
                && (force.abs() < static_force
                    || (force.abs() == static_force && (static_rate >= 0.0 || force == 0.0)));
            if held {
                let release = if static_rate < 0.0 && force != 0.0 {
                    let t = (static_force - force.abs()) / -static_rate;
                    (t >= 0.0 && t <= h_max).then_some(t)
                } else {
                    None
                };
                let h = release.unwrap_or(h_max);
                if h > 0.0 {
                    self.accumulate(
                        &mut result,
                        offset,
                        h,
                        0.0,
                        0.0,
                        force * h,
                        force * h * h * 0.5,
                        0.0,
                    )?;
                    offset += h;
                }
                if release.is_some() {
                    self.push_event(&mut result, offset, DmtSlipEventKind::SlipStarted, 0.0)?;
                    // Compute the following sliding segment from the actual
                    // threshold even when floating-point evaluation lands just
                    // inside its cone. No synthetic nonzero velocity is used.
                    let left = phase_end - offset;
                    if left > 0.0 {
                        let k0 = kk * (area + area_rate * h);
                        self.slide(
                            &mut result,
                            &mut velocity,
                            &mut offset,
                            left,
                            force.signum(),
                            k0,
                            kinetic_rate,
                            force,
                            gate,
                        )?;
                    }
                }
            } else {
                let direction = if velocity != 0.0 {
                    velocity.signum()
                } else {
                    force.signum()
                };
                if direction == 0.0 {
                    return Err(DmtSlipEventError::UnresolvedEvent);
                }
                if velocity == 0.0 {
                    self.push_event(&mut result, offset, DmtSlipEventKind::SlipStarted, 0.0)?;
                }
                self.slide(
                    &mut result,
                    &mut velocity,
                    &mut offset,
                    h_max,
                    direction,
                    kinetic_force,
                    kinetic_rate,
                    force,
                    gate,
                )?;
            }
            if offset == phase_end && boundary == Some(phase_end) && !contact_boundary_done {
                self.push_event(
                    &mut result,
                    offset,
                    if normal_rate > 0.0 {
                        DmtSlipEventKind::ContactEntered
                    } else {
                        DmtSlipEventKind::ContactLeft
                    },
                    0.0,
                )?;
                contact_boundary_done = true;
            }
        }
        result.applied_tangent_work_j =
            finite(applied_tangent_force_n * result.tangent_increment_m)?;
        result.tangent_kinetic_energy_j = finite(0.5 * self.mass_kg * velocity * velocity)?;
        result.total_energy_j =
            finite(result.tangent_kinetic_energy_j + contact.reversible_energy_j)?;
        let kinetic_change = finite(
            0.5 * self.mass_kg
                * (velocity - previous.tangent_velocity_m_s)
                * (velocity + previous.tangent_velocity_m_s),
        )?;
        result.energy_residual_j =
            finite(kinetic_change + result.friction_dissipation_j - result.applied_tangent_work_j)?;
        let energy_tolerance = self.budget.energy_absolute_tolerance_j
            + self.budget.energy_relative_tolerance
                * (kinetic_change.abs()
                    + result.friction_dissipation_j
                    + result.applied_tangent_work_j.abs());
        if result.energy_residual_j.abs() > energy_tolerance {
            return Err(DmtSlipEventError::EnergyBalance {
                residual_j: result.energy_residual_j,
                tolerance_j: energy_tolerance,
            });
        }
        let momentum = self.mass_kg * (velocity - previous.tangent_velocity_m_s);
        result.impulse_residual_n_s =
            finite(momentum - applied_tangent_force_n * dt_s + result.friction_impulse_n_s)?;
        let impulse_tolerance = self.budget.impulse_absolute_tolerance_n_s
            + self.budget.impulse_relative_tolerance
                * (momentum.abs()
                    + (applied_tangent_force_n * dt_s).abs()
                    + result.friction_impulse_n_s.abs());
        if result.impulse_residual_n_s.abs() > impulse_tolerance {
            return Err(DmtSlipEventError::ImpulseBalance {
                residual_n_s: result.impulse_residual_n_s,
                tolerance_n_s: impulse_tolerance,
            });
        }
        result.state = DmtSlipEventState {
            identity: self.identity,
            indentation_m: next_indentation_m,
            tangent_displacement_m: finite(
                previous.tangent_displacement_m + result.tangent_increment_m,
            )?,
            tangent_velocity_m_s: velocity,
            time_s: time,
            steps,
        };
        poll(gate)?;
        Ok(result)
    }

    #[allow(clippy::too_many_arguments)]
    fn slide(
        &self,
        result: &mut DmtSlipEventStep,
        velocity: &mut f64,
        offset: &mut f64,
        h_max: f64,
        direction: f64,
        k0: f64,
        k1: f64,
        force: f64,
        gate: &CancelGate,
    ) -> Result<(), DmtSlipEventError> {
        poll(gate)?;
        let v0 = *velocity;
        let a0 = finite((force - direction * k0) / self.mass_kg)?;
        let a1 = finite(-direction * k1 / self.mass_kg)?;
        let stop = first_zero(v0, a0 * h_max, 0.5 * a1 * h_max * h_max)?;
        let h = stop.map_or(h_max, |u| h_max * u);
        if h <= 0.0 || *offset + h <= *offset {
            return Err(DmtSlipEventError::UnresolvedEvent);
        }
        let value = |t: f64| v0 + t * (a0 + 0.5 * a1 * t);
        let endpoint = finite(value(h))?;
        self.check_velocity(v0)?;
        self.check_velocity(endpoint)?;
        if a1 != 0.0 {
            let extremum = -a0 / a1;
            if extremum > 0.0 && extremum < h {
                self.check_velocity(value(extremum))?;
            }
        }
        // Positive quadrature avoids cancellation between large terms in the
        // exact cubic dissipation polynomial. Its nodes lie strictly inside a
        // sign-preserving segment; an unlocated reversal refuses.
        let t1 = h * 0.211_324_865_405_187_13;
        let t2 = h * 0.788_675_134_594_812_9;
        let w1 = direction * value(t1);
        let w2 = direction * value(t2);
        let c1 = k0 + k1 * t1;
        let c2 = k0 + k1 * t2;
        if [w1, w2, c1, c2].iter().any(|x| !x.is_finite() || *x < 0.0) {
            return Err(DmtSlipEventError::UnresolvedEvent);
        }
        let slip = 0.5 * h * (w1 + w2);
        let impulse = direction * h * (k0 + 0.5 * k1 * h);
        let moment = direction * h * h * (0.5 * k0 + k1 * h / 3.0);
        let loss = 0.5 * h * (c1 * w1 + c2 * w2);
        self.accumulate(
            result,
            *offset,
            h,
            direction * slip,
            slip,
            impulse,
            moment,
            loss,
        )?;
        *offset += h;
        if stop.is_some() {
            let residual = self.mass_kg * endpoint;
            let tolerance = self.budget.impulse_absolute_tolerance_n_s
                + self.budget.impulse_relative_tolerance
                    * (self.mass_kg * v0.abs() + force.abs() * h + impulse.abs());
            if residual.abs() > tolerance {
                return Err(DmtSlipEventError::ImpulseBalance {
                    residual_n_s: residual,
                    tolerance_n_s: tolerance,
                });
            }
            *velocity = 0.0;
            self.push_event(result, *offset, DmtSlipEventKind::Stopped, endpoint)?;
        } else {
            if direction * endpoint < 0.0 {
                return Err(DmtSlipEventError::UnresolvedEvent);
            }
            *velocity = endpoint;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn accumulate(
        &self,
        r: &mut DmtSlipEventStep,
        offset: f64,
        _h: f64,
        dx: f64,
        slip: f64,
        impulse: f64,
        local_moment: f64,
        loss: f64,
    ) -> Result<(), DmtSlipEventError> {
        r.tangent_increment_m = finite(r.tangent_increment_m + dx)?;
        r.accumulated_slip_m = finite(r.accumulated_slip_m + slip)?;
        r.friction_impulse_n_s = finite(r.friction_impulse_n_s + impulse)?;
        r.friction_time_moment_n_s2 =
            finite(r.friction_time_moment_n_s2 + offset * impulse + local_moment)?;
        r.friction_dissipation_j = finite(r.friction_dissipation_j + loss)?;
        r.segments += 1;
        if r.accumulated_slip_m > self.budget.maximum_accumulated_slip_m {
            return Err(DmtSlipEventError::Range("accumulated slip"));
        }
        Ok(())
    }
    fn push_event(
        &self,
        r: &mut DmtSlipEventStep,
        offset: f64,
        kind: DmtSlipEventKind,
        residual: f64,
    ) -> Result<(), DmtSlipEventError> {
        if r.event_count >= self.budget.maximum_events {
            return Err(DmtSlipEventError::EventBudget);
        }
        r.events[r.event_count] = Some(DmtSlipEvent {
            offset_s: offset,
            tangent_increment_m: r.tangent_increment_m,
            kind,
            stop_velocity_residual_m_s: residual,
        });
        r.event_count += 1;
        Ok(())
    }
    fn check_normal(&self, x: f64) -> Result<DmtResponse, DmtSlipEventError> {
        if x < self.budget.minimum_indentation_m {
            return Err(DmtSlipEventError::Range("normal approach"));
        }
        let response = self.normal.evaluate(x)?;
        self.check_force(response.normal_force_n)?;
        self.check_force(self.shear.law().static_shear_strength_pa() * response.contact_area_m2)?;
        self.check_force(self.shear.law().kinetic_shear_strength_pa() * response.contact_area_m2)?;
        Ok(response)
    }
    fn check_force(&self, x: f64) -> Result<(), DmtSlipEventError> {
        if !x.is_finite() || x.abs() > self.budget.maximum_force_n {
            Err(DmtSlipEventError::Range("force"))
        } else {
            Ok(())
        }
    }
    fn check_velocity(&self, v: f64) -> Result<(), DmtSlipEventError> {
        if !v.is_finite() || v.abs() > self.budget.maximum_tangent_speed_m_s {
            Err(DmtSlipEventError::Range("tangent speed"))
        } else {
            Ok(())
        }
    }
}

// First strictly positive root on [0,1]. Scaling prevents overflow in the
// discriminant; the cancellation-resistant q form retains small stop times.
fn first_zero(a: f64, b: f64, c: f64) -> Result<Option<f64>, DmtSlipEventError> {
    if [a, b, c].iter().any(|x| !x.is_finite()) {
        return Err(DmtSlipEventError::UnresolvedEvent);
    }
    let scale = a.abs().max(b.abs()).max(c.abs());
    if scale == 0.0 {
        return Ok(None);
    }
    let (a, b, c) = (a / scale, b / scale, c / scale);
    let admitted = |u: f64| {
        if u.is_finite() && u > 0.0 && u <= 1.0 {
            Some(u)
        } else {
            None
        }
    };
    if a == 0.0 {
        return Ok(if c == 0.0 { None } else { admitted(-b / c) });
    }
    if c == 0.0 {
        return Ok(if b == 0.0 { None } else { admitted(-a / b) });
    }
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        if discriminant.abs() <= 32.0 * f64::EPSILON * (b * b + (4.0 * a * c).abs()) {
            return Err(DmtSlipEventError::UnresolvedEvent);
        }
        return Ok(None);
    }
    let q = -0.5 * (b + det::sqrt(discriminant).copysign(b));
    let u = admitted(q / c);
    let v = if q != 0.0 { admitted(a / q) } else { None };
    Ok(match (u, v) {
        (Some(u), Some(v)) => Some(u.min(v)),
        (u, None) => u,
        (None, v) => v,
    })
}
fn finite(v: f64) -> Result<f64, DmtSlipEventError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(DmtSlipEventError::Range("finite representation"))
    }
}
fn poll(gate: &CancelGate) -> Result<(), DmtSlipEventError> {
    if gate.is_requested() {
        Err(DmtSlipEventError::Cancelled)
    } else {
        Ok(())
    }
}
