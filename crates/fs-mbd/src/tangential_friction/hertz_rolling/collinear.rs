//! Exact force/rate phases with physical slip arrest and rolling-induced restart.
//!
//! The admitted invariant subspace is longitudinal shear plus lateral rolling.
//! Finite-body full-inertia response is retained, but poses, contact points,
//! normal load, material parameters and frequency are frozen for the window.
//! This is a contact-momentum trajectory, not a moving-geometry trajectory.
//! There is no external force or normal impulse in this owner.
//!
//! During slip the original Hertz pure-translation force is `-mu_k*N*sign(u)`.
//! Rolling obeys the original material torque at every instant. An actual zero
//! of shear velocity ends a phase; only there is static *acceleration* balance
//! tested. Thus static capacity cannot arrest a moving contact early. When
//! rolling demands more static force than available, slip starts in the actual
//! outgoing acceleration direction. The native momenta are never projected.

use super::{
    PreparedRigidHertzRollingPair, RigidHertzRollingContact, RigidHertzRollingControls,
    RigidHertzRollingError, finite, rates, sphere_card,
};
use crate::tangential_friction::hertz::multiple::RigidHertzIslandControls;
use crate::tangential_friction::hertz::{
    apply_pair_wrench, points, relative_point, relative_spin, scaled_response,
};
use crate::tangential_friction::rolling::RigidRollingHistory;
use crate::tangential_friction::{
    RigidContactBinding, RigidFrictionError, magnitude, poll, totals,
};
use crate::{DynamicsError, MassProperties, RigidBodyState, Vec3};
use fs_exec::CancelGate;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::viscoelastic::path::{
    ExponentialPathError, ExponentialRatePath, ViscoelasticSpherePathStep,
};
use fs_tribo::sliding_patch::contact::HertzWrenchContact;

#[derive(Clone, Debug, PartialEq)]
struct StoppedCheckpoint {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    body_indices: [usize; 2],
    arms_body: [Vec3; 2],
    frame: TangentFrame,
}

/// Original material history plus an opaque, physically located shear-stop
/// checkpoint. A new raw history never treats a small nonzero velocity as rest.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct RigidHertzRollingCollinearHistory {
    rolling: RigidRollingHistory,
    stopped: Option<StoppedCheckpoint>,
}
impl RigidHertzRollingCollinearHistory {
    #[must_use]
    pub fn new(rolling: RigidRollingHistory) -> Self {
        Self {
            rolling,
            stopped: None,
        }
    }
    #[must_use]
    pub const fn rolling_history(&self) -> &RigidRollingHistory {
        &self.rolling
    }
    #[must_use]
    pub const fn is_stopped(&self) -> bool {
        self.stopped.is_some()
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RigidHertzRollingCollinearControls {
    pub mechanical: RigidHertzRollingControls,
    /// At most4096 actual slip/stick phases in one atomic window.
    pub maximum_phases: usize,
    /// At most256 bracketed root iterations per stop event.
    pub maximum_root_iterations: usize,
    /// Absolute bracket width for a numerically located zero [s]. Its native
    /// velocity residual must also meet the original mechanical tolerance.
    pub time_tolerance_s: f64,
}
impl RigidHertzRollingCollinearControls {
    pub fn validate(self) -> Result<(), RigidHertzRollingCollinearError> {
        RigidHertzIslandControls {
            mechanical: self.mechanical.mechanical,
            maximum_contacts: 1,
            maximum_matrix_terms: 9,
            maximum_sweeps: 1,
            maximum_evaluations: 1,
        }
        .validate()
        .map_err(|_| invalid("invalid original mechanical controls"))?;
        if !(1..=4096).contains(&self.maximum_phases)
            || !(1..=256).contains(&self.maximum_root_iterations)
            || !self.time_tolerance_s.is_finite()
            || self.time_tolerance_s <= 0.0
            || [
                self.mechanical.maximum_rolling_moment_nm,
                self.mechanical.rolling_moment_tolerance_nm,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(invalid(
                "bounded phase/root work and finite physical controls required",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RigidHertzRollingCollinearMode {
    SlidingPositive,
    SlidingNegative,
    Sticking,
}

/// One continuous original-law phase. Rolling moment in the original receipt
/// is its endpoint value; the separate impulse integrates the entire path.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRollingCollinearPhase {
    pub start_time_s: f64,
    pub duration_s: f64,
    pub mode: RigidHertzRollingCollinearMode,
    pub reached_stop: bool,
    pub root_iterations: usize,
    pub root_bracket_s: f64,
    pub shear_velocity_before_m_s: f64,
    pub shear_velocity_after_m_s: f64,
    pub rolling_rate_before_rad_s: f64,
    pub rolling_rate_after_rad_s: f64,
    pub shear_displacement_m: f64,
    pub rolling_angle_rad: f64,
    pub shear_impulse_ns: f64,
    pub rolling_moment_impulse_nms: f64,
    pub maximum_shear_force_n: f64,
    pub surface_heat_j: f64,
    pub bulk_heat_j: f64,
    pub native_shear_velocity_after_m_s: f64,
    pub native_rolling_rate_after_rad_s: f64,
    pub maximum_velocity_residual_m_s: f64,
    pub kinetic_energy_change_j: f64,
    pub energy_residual_j: f64,
    pub rolling: ViscoelasticSpherePathStep,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRollingCollinearReport {
    pub phases: Vec<RigidHertzRollingCollinearPhase>,
    pub deforming_body: usize,
    pub rolling_temperature_k: f64,
    pub rolling_excitation_frequency_hz: f64,
    pub duration_s: f64,
    pub scaled_mobility_per_kg: [[f64; 2]; 2],
    pub surface_heat_j: f64,
    pub bulk_heat_j: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub total_energy_residual_j: f64,
    pub history_heat_residual_j: f64,
    pub force_impulse_world_ns: Vec3,
    pub rolling_angular_impulse_world_nms: Vec3,
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub maximum_velocity_residual_m_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRollingCollinearCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub history: RigidHertzRollingCollinearHistory,
    pub report: RigidHertzRollingCollinearReport,
}

/// Immutable geometry, original inertia/material and the accepted stop history.
#[derive(Clone, Debug)]
pub struct PreparedRigidHertzRollingCollinear<'a> {
    pair: PreparedRigidHertzRollingPair<'a>,
    history: RigidHertzRollingCollinearHistory,
    controls: RigidHertzRollingCollinearControls,
    response: [[f64; 2]; 2],
    gamma: f64,
}

#[derive(Debug)]
pub enum RigidHertzRollingCollinearError {
    Invalid(&'static str),
    Pair(RigidHertzRollingError),
    Path(ExponentialPathError),
    PhaseBudget,
    RootBudget,
    NoProgress,
    EffortBudget {
        required_force_n: f64,
        required_rolling_moment_nm: f64,
    },
    PhaseAudit(Box<RigidHertzRollingCollinearPhase>),
    Audit(Box<RigidHertzRollingCollinearReport>),
}
impl From<RigidHertzRollingError> for RigidHertzRollingCollinearError {
    fn from(e: RigidHertzRollingError) -> Self {
        Self::Pair(e)
    }
}
impl From<ExponentialPathError> for RigidHertzRollingCollinearError {
    fn from(e: ExponentialPathError) -> Self {
        Self::Path(e)
    }
}
impl From<RigidFrictionError> for RigidHertzRollingCollinearError {
    fn from(e: RigidFrictionError) -> Self {
        Self::Pair(e.into())
    }
}
impl From<DynamicsError> for RigidHertzRollingCollinearError {
    fn from(e: DynamicsError) -> Self {
        Self::Pair(e.into())
    }
}
impl From<crate::tangential_friction::hertz::RigidHertzError> for RigidHertzRollingCollinearError {
    fn from(e: crate::tangential_friction::hertz::RigidHertzError) -> Self {
        Self::Pair(e.into())
    }
}
impl core::fmt::Display for RigidHertzRollingCollinearError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "collinear Hertz/rolling trajectory: {self:?}")
    }
}
impl std::error::Error for RigidHertzRollingCollinearError {}
fn invalid(s: &'static str) -> RigidHertzRollingCollinearError {
    RigidHertzRollingCollinearError::Invalid(s)
}

impl<'a> PreparedRigidHertzRollingCollinear<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        history: &RigidHertzRollingCollinearHistory,
        contact: RigidHertzRollingContact<'a>,
        dt_s: f64,
        controls: RigidHertzRollingCollinearControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidHertzRollingCollinearError> {
        poll(gate)?;
        controls.validate()?;
        let pair = PreparedRigidHertzRollingPair::new(
            bodies,
            properties,
            &history.rolling,
            contact,
            dt_s,
            controls.mechanical,
            gate,
        )?;
        let response = [
            [pair.mobility[0][0], pair.mobility[0][4]],
            [pair.mobility[4][0], pair.mobility[4][4]],
        ];
        if response[0][0] <= 0.0
            || response[1][1] <= 0.0
            || finite(response[1][1] - response[0][1] * (response[1][0] / response[0][0]))? <= 0.0
        {
            return Err(invalid("positive original shear/rolling response required"));
        }
        let initial_point = relative_point(pair.core.points_before);
        let initial_spin = relative_spin(pair.core.points_before);
        let point_roundoff = 128.0 * f64::EPSILON * magnitude(initial_point);
        let spin_roundoff = 128.0 * f64::EPSILON * magnitude(initial_spin);
        if initial_point.dot(pair.core.axes[1]).abs() > point_roundoff
            || initial_point.dot(pair.core.axes[2]).abs() > point_roundoff
            || initial_spin.dot(pair.core.axes[0]).abs() > spin_roundoff
            || initial_spin.dot(pair.core.axes[2]).abs() > spin_roundoff
        {
            return Err(invalid(
                "initial motion leaves the invariant longitudinal/lateral-rolling subspace",
            ));
        }
        for i in [1, 2, 3] {
            for j in [0, 4] {
                let bound =
                    128.0 * f64::EPSILON * pair.mobility[i][i].sqrt() * pair.mobility[j][j].sqrt();
                if pair.mobility[i][j].abs() > bound {
                    return Err(invalid(
                        "original full inertia couples shear/rolling to an excluded contact port",
                    ));
                }
            }
        }
        let normal = scaled_response(
            pair.core.before,
            pair.core.pair_properties,
            pair.core.arms,
            [pair.core.axes[2], pair.core.axes[0], pair.core.axes[1]],
            [0.0, 0.0, pair.radius_m],
        )?;
        for j in [1, 2] {
            let bound = 128.0 * f64::EPSILON * normal[0][0].sqrt() * normal[j][j].sqrt();
            if normal[0][j].abs() > bound {
                return Err(invalid(
                    "selected friction ports induce unresolved normal motion",
                ));
            }
        }
        if let Some(stop) = &history.stopped {
            if stop.bodies != pair.core.before
                || stop.properties != pair.core.pair_properties
                || stop.body_indices != [contact.binding.body_a, contact.binding.body_b]
                || stop.arms_body != [contact.binding.arm_a_body, contact.binding.arm_b_body]
                || stop.frame != contact.frame
                || pair.rate_before[0].abs() > controls.mechanical.mechanical.velocity_tolerance_m_s
            {
                return Err(invalid(
                    "stopped continuation requires its exact native bodies and contact geometry",
                ));
            }
        }
        let gamma = finite(
            (sphere_card(contact.rolling)?.rolling_time_s() * pair.normal_force_n())
                / pair.radius_m,
        )?;
        if gamma == 0.0
            && sphere_card(contact.rolling)?.rolling_time_s() > 0.0
            && pair.normal_force_n() > 0.0
        {
            return Err(invalid("original material rolling damping underflow"));
        }
        poll(gate)?;
        Ok(Self {
            pair,
            history: history.clone(),
            controls,
            response,
            gamma,
        })
    }
    #[must_use]
    pub const fn binding(&self) -> RigidContactBinding {
        self.pair.contact.binding
    }
    #[must_use]
    pub fn duration_s(&self) -> f64 {
        self.pair.duration_s()
    }
    #[must_use]
    pub fn contact_points_world(&self) -> [Vec3; 2] {
        self.pair.core.contact_points_world()
    }
    #[must_use]
    pub const fn sphere_radius_m(&self) -> f64 {
        self.pair.radius_m
    }
    #[must_use]
    pub const fn deforming_body(&self) -> usize {
        self.pair.deforming_body
    }
    #[must_use]
    pub fn normal_force_n(&self) -> f64 {
        self.pair.normal_force_n()
    }
    #[must_use]
    pub fn contact_radius_m(&self) -> f64 {
        self.pair.contact_radius_m()
    }
    pub fn advance(
        &self,
        law: &HertzWrenchContact,
        gate: &CancelGate,
    ) -> Result<RigidHertzRollingCollinearCandidate, RigidHertzRollingCollinearError> {
        self.advance_at_temperature(law, self.pair.contact.temperature_k, gate)
    }
    #[allow(clippy::too_many_lines)]
    pub fn advance_at_temperature(
        &self,
        law: &HertzWrenchContact,
        temperature_k: f64,
        gate: &CancelGate,
    ) -> Result<RigidHertzRollingCollinearCandidate, RigidHertzRollingCollinearError> {
        poll(gate)?;
        let parameters = law.parameters();
        if parameters.normal_force_n.to_bits() != self.normal_force_n().to_bits()
            || parameters.contact_radius_m.to_bits() != self.contact_radius_m().to_bits()
        {
            return Err(invalid(
                "prepared original normal load and radius must remain unchanged",
            ));
        }
        let static_capacity = finite(parameters.static_mu * parameters.normal_force_n)?;
        let kinetic_capacity = finite(parameters.kinetic_mu * parameters.normal_force_n)?;
        let [[a, b], [_b, c]] = self.response;
        let dt = self.duration_s();
        let radius = self.pair.radius_m;
        let mut current = self.pair.core.before;
        let mut rolling_history = self.history.rolling.clone();
        let mut stopped = self.history.is_stopped() || self.pair.rate_before[0] == 0.0;
        let mut elapsed = 0.0;
        let mut phases = Vec::new();
        let mut total_shear_impulse = 0.0;
        let mut total_rolling_impulse = 0.0;
        let m = self.controls.mechanical.mechanical;
        while elapsed < dt {
            poll(gate)?;
            if phases.len() == self.controls.maximum_phases {
                return Err(RigidHertzRollingCollinearError::PhaseBudget);
            }
            let remaining = dt - elapsed;
            let p = points(current, self.pair.core.pair_properties, self.pair.core.arms)?;
            let native_rates = rates(p, self.pair.core.axes, self.contact_radius_m(), radius);
            let u = if stopped { 0.0 } else { native_rates[0] };
            let r = native_rates[4];
            let required_static = finite((b / a) * self.gamma * r)?;
            let sticking = stopped && required_static.abs() <= static_capacity;
            let sign = if u != 0.0 {
                u.signum()
            } else {
                (-b * self.gamma * r).signum()
            };
            let force = if sticking {
                0.0
            } else {
                -kinetic_capacity * sign
            };
            let drive = if sticking { 0.0 } else { b * force };
            let decay = if sticking {
                self.gamma * (c - b * (b / a))
            } else {
                self.gamma * c
            };
            let mode = if sticking {
                RigidHertzRollingCollinearMode::Sticking
            } else if sign > 0.0 {
                RigidHertzRollingCollinearMode::SlidingPositive
            } else {
                RigidHertzRollingCollinearMode::SlidingNegative
            };
            let stop = if sticking {
                None
            } else {
                self.first_stop(u, r, force, sign, remaining, gate)?
            };
            let duration = stop.as_ref().map_or(remaining, |root| root.time);
            if !duration.is_finite() || duration <= 0.0 || elapsed + duration <= elapsed {
                return Err(RigidHertzRollingCollinearError::NoProgress);
            }
            let path = ExponentialRatePath::new(r, drive, decay, duration, || gate.is_requested())?;
            let angular_path =
                ExponentialRatePath::new(r / radius, drive / radius, decay, duration, || {
                    gate.is_requested()
                })?;
            let rolling = self.pair.contact.rolling.advance_exponential_sphere_path(
                self.pair.contact.patch,
                self.pair.contact.interface,
                &angular_path,
                0.0,
                temperature_k,
                self.pair.contact.excitation_frequency_hz,
                self.pair.contact.work,
                rolling_history.state(),
                || gate.is_requested(),
            )?;
            let (shear_impulse, travel, shear_after, maximum_force) = if sticking {
                (
                    finite((b / a) * self.gamma * path.rate_integral())?,
                    0.0,
                    0.0,
                    required_static.abs(),
                )
            } else {
                (
                    finite(force * duration)?,
                    finite((-b * self.gamma).mul_add(
                        path.integrated_displacement(),
                        u * duration + 0.5 * a * force * duration * duration,
                    ))?,
                    finite(
                        (-b * self.gamma).mul_add(path.rate_integral(), u + a * force * duration),
                    )?,
                    force.abs(),
                )
            };
            let surface_heat = if sticking {
                0.0
            } else {
                finite(-force * travel)?
            };
            if surface_heat < 0.0 {
                return Err(invalid(
                    "resolved kinetic phase has negative original surface dissipation",
                ));
            }
            let maximum_moment =
                finite(self.gamma * radius * r.abs().max(path.endpoint_rate().abs()))?;
            if maximum_force > m.contact.maximum_force_n
                || maximum_moment > self.controls.mechanical.maximum_rolling_moment_nm
            {
                return Err(RigidHertzRollingCollinearError::EffortBudget {
                    required_force_n: maximum_force,
                    required_rolling_moment_nm: maximum_moment,
                });
            }
            let impulse = self.pair.core.axes[0].scale(shear_impulse);
            let couple = self.pair.core.axes[1].scale(rolling.rolling_moment_impulse_nms);
            let after = apply_pair_wrench(
                current,
                self.pair.core.pair_properties,
                self.pair.core.arms,
                impulse,
                couple,
                gate,
            )?;
            let after_points = points(after, self.pair.core.pair_properties, self.pair.core.arms)?;
            let after_rates = rates(
                after_points,
                self.pair.core.axes,
                self.contact_radius_m(),
                radius,
            );
            let excluded_spin = relative_spin(after_points)
                .dot(self.pair.core.axes[0])
                .abs()
                .max(
                    relative_spin(after_points)
                        .dot(self.pair.core.axes[2])
                        .abs(),
                );
            let velocity_residual = (after_rates[0] - shear_after)
                .abs()
                .max((after_rates[4] - path.endpoint_rate()).abs())
                .max(after_rates[1].abs())
                .max(after_rates[2].abs())
                .max(after_rates[3].abs())
                .max(
                    relative_point(after_points)
                        .dot(self.pair.core.axes[2])
                        .abs(),
                );
            let energy_before = totals(&current, &self.pair.core.pair_properties, gate)?.0;
            let energy_after = totals(&after, &self.pair.core.pair_properties, gate)?.0;
            let bulk_heat = rolling.rolling.dissipation.total_heat_j;
            let phase = RigidHertzRollingCollinearPhase {
                start_time_s: elapsed,
                duration_s: duration,
                mode,
                reached_stop: stop.is_some(),
                root_iterations: stop.as_ref().map_or(0, |s| s.iterations),
                root_bracket_s: stop.as_ref().map_or(0.0, |s| s.bracket),
                shear_velocity_before_m_s: u,
                shear_velocity_after_m_s: shear_after,
                rolling_rate_before_rad_s: r / radius,
                rolling_rate_after_rad_s: path.endpoint_rate() / radius,
                shear_displacement_m: travel,
                rolling_angle_rad: angular_path.rate_integral(),
                shear_impulse_ns: shear_impulse,
                rolling_moment_impulse_nms: rolling.rolling_moment_impulse_nms,
                maximum_shear_force_n: maximum_force,
                surface_heat_j: surface_heat,
                bulk_heat_j: bulk_heat,
                native_shear_velocity_after_m_s: after_rates[0],
                native_rolling_rate_after_rad_s: after_rates[4] / radius,
                maximum_velocity_residual_m_s: velocity_residual,
                kinetic_energy_change_j: energy_after - energy_before,
                energy_residual_j: energy_after - energy_before + surface_heat + bulk_heat,
                rolling,
            };
            if !phase.energy_residual_j.is_finite()
                || phase.energy_residual_j.abs() > m.energy_tolerance_j
                || !velocity_residual.is_finite()
                || velocity_residual > m.velocity_tolerance_m_s
                || (after_rates[4] - path.endpoint_rate()).abs() / radius
                    > m.angular_velocity_tolerance_rad_s
                || excluded_spin > m.angular_velocity_tolerance_rad_s
                || (phase.reached_stop && after_rates[0].abs() > m.velocity_tolerance_m_s)
            {
                return Err(RigidHertzRollingCollinearError::PhaseAudit(Box::new(phase)));
            }
            rolling_history = RigidRollingHistory::from_native_step(&phase.rolling.rolling);
            stopped = sticking || phase.reached_stop;
            total_shear_impulse = finite(total_shear_impulse + shear_impulse)?;
            total_rolling_impulse =
                finite(total_rolling_impulse + phase.rolling_moment_impulse_nms)?;
            elapsed = if duration == remaining {
                dt
            } else {
                elapsed + duration
            };
            current = after;
            phases.push(phase);
        }
        // A final single native summed impulse is the published body endpoint;
        // intermediate states exist only for path/domain/event verification.
        let impulse = self.pair.core.axes[0].scale(total_shear_impulse);
        let couple = self.pair.core.axes[1].scale(total_rolling_impulse);
        let after = apply_pair_wrench(
            self.pair.core.before,
            self.pair.core.pair_properties,
            self.pair.core.arms,
            impulse,
            couple,
            gate,
        )?;
        let before_totals = self.pair.core.totals_before;
        let after_totals = totals(&after, &self.pair.core.pair_properties, gate)?;
        let final_points = points(after, self.pair.core.pair_properties, self.pair.core.arms)?;
        let final_rates = rates(
            final_points,
            self.pair.core.axes,
            self.contact_radius_m(),
            radius,
        );
        let final_excluded_spin = relative_spin(final_points)
            .dot(self.pair.core.axes[0])
            .abs()
            .max(
                relative_spin(final_points)
                    .dot(self.pair.core.axes[2])
                    .abs(),
            );
        let last = phases
            .last()
            .expect("positive admitted duration produced a phase");
        let final_velocity_residual = (final_rates[0] - last.shear_velocity_after_m_s)
            .abs()
            .max((final_rates[4] - radius * last.rolling_rate_after_rad_s).abs())
            .max(final_rates[1].abs())
            .max(final_rates[2].abs())
            .max(final_rates[3].abs())
            .max(
                relative_point(final_points)
                    .dot(self.pair.core.axes[2])
                    .abs(),
            );
        if !final_velocity_residual.is_finite()
            || final_velocity_residual > m.velocity_tolerance_m_s
            || (stopped && final_rates[0].abs() > m.velocity_tolerance_m_s)
            || (final_rates[4] / radius - last.rolling_rate_after_rad_s).abs()
                > m.angular_velocity_tolerance_rad_s
            || final_excluded_spin > m.angular_velocity_tolerance_rad_s
        {
            return Err(invalid(
                "summed native endpoint does not close the resolved contact trajectory",
            ));
        }
        let surface_heat = phases.iter().map(|p| p.surface_heat_j).sum::<f64>();
        let bulk_heat = phases.iter().map(|p| p.bulk_heat_j).sum::<f64>();
        let maximum_residual = phases
            .iter()
            .map(|p| p.maximum_velocity_residual_m_s)
            .fold(final_velocity_residual, f64::max);
        let report = RigidHertzRollingCollinearReport {
            phases,
            deforming_body: self.deforming_body(),
            rolling_temperature_k: temperature_k,
            rolling_excitation_frequency_hz: self.pair.contact.excitation_frequency_hz,
            duration_s: dt,
            scaled_mobility_per_kg: self.response,
            surface_heat_j: surface_heat,
            bulk_heat_j: bulk_heat,
            kinetic_before_j: before_totals.0,
            kinetic_after_j: after_totals.0,
            total_energy_residual_j: after_totals.0 - before_totals.0 + surface_heat + bulk_heat,
            history_heat_residual_j: rolling_history.state().cumulative_heat_j()
                - self.history.rolling.state().cumulative_heat_j()
                - bulk_heat,
            force_impulse_world_ns: impulse,
            rolling_angular_impulse_world_nms: couple,
            linear_momentum_change_ns: after_totals.1.sub(before_totals.1),
            angular_momentum_change_nms: after_totals.2.sub(before_totals.2),
            maximum_velocity_residual_m_s: maximum_residual,
        };
        if [
            report.total_energy_residual_j,
            report.history_heat_residual_j,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v.abs() > m.energy_tolerance_j)
            || magnitude(report.linear_momentum_change_ns) > m.linear_momentum_tolerance_ns
            || magnitude(report.angular_momentum_change_nms) > m.angular_momentum_tolerance_nms
        {
            return Err(RigidHertzRollingCollinearError::Audit(Box::new(report)));
        }
        let history = RigidHertzRollingCollinearHistory {
            rolling: rolling_history,
            stopped: stopped.then_some(StoppedCheckpoint {
                bodies: after,
                properties: self.pair.core.pair_properties,
                body_indices: [self.binding().body_a, self.binding().body_b],
                arms_body: [self.binding().arm_a_body, self.binding().arm_b_body],
                frame: self.pair.contact.frame,
            }),
        };
        poll(gate)?;
        Ok(RigidHertzRollingCollinearCandidate {
            bodies_after: after,
            history,
            report,
        })
    }

    fn slip_at(
        &self,
        u: f64,
        r: f64,
        force: f64,
        t: f64,
        gate: &CancelGate,
    ) -> Result<(f64, f64, f64), RigidHertzRollingCollinearError> {
        let [[a, b], [_b, c]] = self.response;
        let path =
            ExponentialRatePath::new(r, b * force, self.gamma * c, t, || gate.is_requested())?;
        let arithmetic_scale =
            u.abs() + (a * force * t).abs() + (b * self.gamma * path.rate_integral()).abs();
        let root_tolerance = self
            .controls
            .mechanical
            .mechanical
            .velocity_tolerance_m_s
            .min(128.0 * f64::EPSILON * arithmetic_scale);
        Ok((
            finite((-b * self.gamma).mul_add(path.rate_integral(), u + a * force * t))?,
            finite(a * force - b * self.gamma * path.endpoint_rate())?,
            root_tolerance,
        ))
    }

    /// For constant Coulomb force, du/dt is a constant plus one exponential.
    /// Its limiting sign opposes slip, so there is at most one positive stop
    /// after any initial maximum. The bracket excludes the departure zero.
    fn first_stop(
        &self,
        u: f64,
        r: f64,
        force: f64,
        sign: f64,
        horizon: f64,
        gate: &CancelGate,
    ) -> Result<Option<StopRoot>, RigidHertzRollingCollinearError> {
        let (end, _, _) = self.slip_at(u, r, force, horizon, gate)?;
        if sign * end > 0.0 || (u == 0.0 && force == 0.0) {
            return Ok(None);
        }
        if end == 0.0 {
            return Ok(Some(StopRoot {
                time: horizon,
                iterations: 0,
                bracket: 0.0,
            }));
        }
        let [[a, b], [_b, c]] = self.response;
        let mut lo = 0.0;
        if u == 0.0 {
            let initial = sign * (a * force - b * self.gamma * r);
            let limiting = sign * (a - b * (b / c)) * force;
            if initial <= 0.0 || limiting >= 0.0 || self.gamma == 0.0 {
                return Err(invalid(
                    "unresolved departure from an inadmissible static state",
                ));
            }
            let ratio = initial / (-limiting);
            let logarithm = if ratio.is_finite() {
                ratio.ln_1p()
            } else {
                initial.ln() - (-limiting).ln()
            };
            lo = finite(logarithm / (self.gamma * c))?;
            if lo <= 0.0 || lo >= horizon {
                return Err(invalid("unresolved outgoing slip extremum"));
            }
        }
        let mut hi = horizon;
        let initial_derivative = a * force - b * self.gamma * r;
        let linear_guess = -u / initial_derivative;
        let mut t = if linear_guess.is_finite() && linear_guess > lo && linear_guess < hi {
            linear_guess
        } else {
            lo + 0.5 * (hi - lo)
        };
        for iteration in 1..=self.controls.maximum_root_iterations {
            poll(gate)?;
            let (value, derivative, root_tolerance) = self.slip_at(u, r, force, t, gate)?;
            if value == 0.0
                || (hi - lo <= self.controls.time_tolerance_s && value.abs() <= root_tolerance)
            {
                return Ok(Some(StopRoot {
                    time: t,
                    iterations: iteration,
                    bracket: hi - lo,
                }));
            }
            if sign * value > 0.0 {
                lo = t;
            } else {
                hi = t;
            }
            let newton = t - value / derivative;
            let next = if newton.is_finite() && newton > lo && newton < hi {
                newton
            } else {
                lo + 0.5 * (hi - lo)
            };
            if next == lo || next == hi || next == t {
                if value.abs() <= root_tolerance && hi - lo <= self.controls.time_tolerance_s {
                    return Ok(Some(StopRoot {
                        time: t,
                        iterations: iteration,
                        bracket: hi - lo,
                    }));
                }
                return Err(RigidHertzRollingCollinearError::NoProgress);
            }
            t = next;
        }
        Err(RigidHertzRollingCollinearError::RootBudget)
    }
}

struct StopRoot {
    time: f64,
    iterations: usize,
    bracket: f64,
}
