//! Translating fiber supports, their inertia and the moving material crossing.
//!
//! In the support frame `mu u_tt + elastic/loss terms = contact - mu X_tt`.
//! The original Galerkin mass projection supplies the inertial load; no hair
//! mass, spring or damping is invented here. A fixed laboratory crossing samples
//! `s = s_initial - X(t)`. Its material velocity is `X_t + phi(s)^T q_t`, NOT
//! the time derivative of `phi(s(t))^T q` (which follows the moving crossing).
//!
//! Straight parallel fibers and prescribed normal loads remain this rung's
//! assumptions. Normal patches and rosin histories are not transported. The
//! native mechanical ledger is a SUPPORT-FRAME balance: external work includes
//! the inertial load. `FiberSupportFrame` additionally accounts for rigid
//! translation and momentum cross terms in laboratory energy and support work.
//! Refine time, span modes and width.
use super::super::stations::BowStationMotionLimits;
use super::super::*;
use crate::fiber_ribbon::FiberDirection;
use crate::render::RenderError;
use crate::render::schedule::force::coupled::contact::friction::multiple::surfaces::geometry::{
    FrictionEndpointGeometry, FrictionGeometry, FrictionGeometryConfig,
};

/// Physical stroke and numerical map-update admission; no automatic clamping.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberSupportMotionLimits {
    /// Absolute axial support speed at both step endpoints [m/s].
    pub maximum_speed_m_s: f64,
    /// Absolute acceleration over one mechanical interval [m/s^2].
    pub maximum_acceleration_m_s2: f64,
    /// Conservative travel bound `dt*max(abs(v0),abs(v1))` [m].
    pub maximum_step_m: f64,
    /// Work budget for the original complete contact-map response rebuild.
    pub max_update_terms: usize,
}
impl FiberSupportMotionLimits {
    pub fn validate(self) -> Result<(), ModalFrictionError> {
        BowStationMotionLimits {
            maximum_speed_m_s: self.maximum_speed_m_s,
            maximum_step_m: self.maximum_step_m,
            max_update_terms: self.max_update_terms,
        }
        .validate()?;
        if !self.maximum_acceleration_m_s2.is_finite() || self.maximum_acceleration_m_s2 <= 0.0 {
            return Err(error(
                "fiber support requires a positive finite acceleration limit",
            ));
        }
        Ok(())
    }
}
/// Last accepted support endpoint relative to the original material description.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberSupportKinematics {
    pub velocity_m_s: f64,
    pub displacement_m: f64,
    pub station_fraction: f64,
}
/// Laboratory-frame balance of one accepted step of the original fiber owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberSupportFrame {
    pub sample: u64,
    pub velocity_m_s: f64,
    pub translation_m: f64,
    pub acceleration_m_s2: f64,
    /// `-a sum Gamma_n delta_q_n`, already included in native external work [J].
    pub inertial_coordinate_work_j: f64,
    /// `V sum Gamma_n qdot_n + M V^2/2`, added to native storage [J].
    pub laboratory_kinetic_correction_j: f64,
    /// Native drive work + inertial work + change of kinetic correction [J].
    /// Other applied forces and moving normal supports keep their own accounts.
    pub support_work_j: f64,
    /// Complete accepted mechanical storage plus the kinetic correction [J].
    pub laboratory_mechanical_energy_j: f64,
}
pub(crate) struct FiberSupportState {
    limits: FiberSupportMotionLimits,
    accepted: FiberSupportKinematics,
    participation: Vec<f64>,
    force_start: usize,
    mass_kg: f64,
    last_frame: Option<FiberSupportFrame>,
}
pub(crate) struct FiberSupportStep {
    accepted: FiberSupportKinematics,
    component_start: usize,
    shape: Vec<f64>,
    start_shape: Vec<f64>,
    end_shape: Vec<f64>,
    start_speed_m_s: f64,
    acceleration_m_s2: f64,
    previous_displacement_moment: f64,
    previous_kinetic_correction_j: f64,
    pub(crate) forces: Vec<f64>,
    pub(crate) drives: Vec<FrictionDrive>,
}
impl FiberSupportStep {
    pub(crate) fn endpoint_geometry(
        &self,
        cell: usize,
        left_start: ModalAttachment,
        left_end: ModalAttachment,
    ) -> FrictionEndpointGeometry {
        FrictionEndpointGeometry {
            start: FrictionGeometry {
                left: left_start,
                right: Some(ModalAttachment {
                    component: self.component_start + cell,
                    shapes: self.start_shape.clone(),
                }),
            },
            end: FrictionGeometry {
                left: left_end,
                right: Some(ModalAttachment {
                    component: self.component_start + cell,
                    shapes: self.end_shape.clone(),
                }),
            },
            start_speed_m_s: self.start_speed_m_s,
            end_speed_m_s: self.accepted.velocity_m_s,
        }
    }
    pub(crate) fn receiver(&self, cell: usize) -> ModalAttachment {
        ModalAttachment {
            component: self.component_start + cell,
            shapes: self.shape.clone(),
        }
    }
}
fn error(what: &'static str) -> ModalFrictionError {
    ModalCouplingError::Invalid { what }.into()
}
fn admit_sources(
    source: &MultiBowStringBankSource,
) -> Result<&ribbon::MaterialFiberRibbon, ModalFrictionError> {
    let fiber = source
        .fiber_ribbon()
        .ok_or_else(|| error("support motion requires material fibers"))?;
    let start = fiber.component_start();
    let end = start + fiber.model().resolution().cells;
    // Bow-position motion on the string may carry a normal patch, but travel
    // along accelerating hair still lacks that normal footprint's transport.
    // Keep this restriction local to fiber motion rather than the shared
    // bow-station admission, which also serves fixed material receivers.
    if !source.bow_contact_loads().is_empty()
        || source.stops().iter().any(|stop| {
            stop.body.as_ref().is_some_and(|body| {
                (start..end).contains(&(source.bank().strings().len() + body.component))
            })
        })
    {
        return Err(error(
            "fiber support motion requires transported normal hair contacts",
        ));
    }
    Ok(fiber)
}
impl MultiBowStringBank {
    /// Opt into finite fiber travel and acceleration before playback.
    /// All cells share ONE support. Drive speeds then denote ENDPOINT velocities;
    /// constant acceleration connects accepted endpoints. Friction receives the
    /// interval mean; each fiber receives its mass-projected inertial load and
    /// a moving contact in the original basis, without restarting its waves.
    ///
    /// Bind after source contacts. Only prescribed-load dry friction is admitted;
    /// loaded normals and shear/thermal history transport remain unsupported.
    /// Initial speed declares an already co-moving state, not an impulse. A
    /// velocity step in a low-rate schedule remains a finite-interval acceleration
    /// and may fail the stated acceleration bound. Refine control timing too.
    pub fn with_fiber_support_motion(
        mut self,
        initial_velocity_m_s: f64,
        limits: FiberSupportMotionLimits,
        gate: &CancelGate,
    ) -> Result<Self, RenderError> {
        if gate.is_requested() {
            return Err(RenderError::Coupled(ModalCouplingError::Cancelled));
        }
        limits.validate().map_err(RenderError::Friction)?;
        if self.system.samples_rendered() != 0
            || !initial_velocity_m_s.is_finite()
            || initial_velocity_m_s.abs() > limits.maximum_speed_m_s
            || self
                .bow_stations
                .as_ref()
                .is_some_and(|s| s.fiber_support.is_some())
        {
            return Err(RenderError::Friction(error(
                "fiber motion binds once before playback with an admitted initial speed",
            )));
        }
        self.admit_bow_station_sources()
            .map_err(RenderError::Friction)?;
        let fiber = admit_sources(&self.source).map_err(RenderError::Friction)?;
        let station_fraction = fiber.station_fraction();
        let participation = fiber.model().axial_translation_participation();
        let mass_kg = fiber.model().mass_kg();
        let component_start = fiber.component_start();
        let component_end = component_start
            .checked_add(fiber.model().resolution().cells)
            .ok_or_else(|| RenderError::Friction(error("fiber component span overflow")))?;
        let components = self
            .system
            .components()
            .get(component_start..component_end)
            .ok_or_else(|| RenderError::Friction(error("original fiber components are absent")))?;
        // The retained cell layout contains axial THEN normal modes. Reuse that
        // complete source projection; never infer offsets from axial count alone.
        if participation.iter().any(|v| !v.is_finite())
            || components
                .iter()
                .any(|c| c.modes().len() != participation.len())
        {
            return Err(RenderError::Friction(error(
                "fiber participation does not match the retained components",
            )));
        }
        // Every candidate enforces these original component state limits.
        // Admit the whole frame conversion BEFORE any step, leaving generous
        // floating-point headroom so accepted mechanics cannot fail afterward
        // merely while publishing laboratory work/energy.
        let mut displacement_bound = 0.0;
        let mut momentum_bound = 0.0;
        for component in components {
            let budget = component.budget();
            for gamma in &participation {
                displacement_bound += gamma.abs() * budget.maximum_abs_displacement_m_sqrt_kg;
                momentum_bound += gamma.abs() * budget.maximum_abs_velocity_m_sqrt_kg_per_s;
            }
        }
        let correction_bound = limits.maximum_speed_m_s * momentum_bound
            + 0.5 * mass_kg * limits.maximum_speed_m_s * limits.maximum_speed_m_s;
        let inertial_work_bound = limits.maximum_acceleration_m_s2 * (2.0 * displacement_bound);
        let drive_bound = self
            .source
            .bows()
            .iter()
            .map(|bow| {
                bow.friction.maximum_force_n
                    * limits.maximum_speed_m_s
                    * self.system.sample_period_s()
            })
            .sum::<f64>();
        // The joint owner bounds total storage, including admitted obstacles.
        // Reserve each separately reported storage field too, independently
        // of its partition into network, normal, residual and shear accounts.
        let account_bound = 4.0 * self.source.bank().options().network.maximum_total_energy_j
            + drive_bound
            + inertial_work_bound
            + 2.0 * correction_bound;
        if !account_bound.is_finite() || account_bound > f64::MAX / 16.0 {
            return Err(RenderError::Friction(error(
                "fiber laboratory frame exceeds representable state/work bounds",
            )));
        }
        let force_start = self.system.components()[..fiber.component_start()]
            .iter()
            .map(|c| c.modes().len())
            .sum();
        if self.bow_stations.is_none() {
            self = self.with_bow_station_motion(
                BowStationMotionLimits {
                    maximum_speed_m_s: limits.maximum_speed_m_s,
                    maximum_step_m: limits.maximum_step_m,
                    max_update_terms: limits.max_update_terms,
                },
                gate,
            )?;
        }
        let state = self.bow_stations.as_mut().expect("admitted geometry");
        let max_update_terms = state.limits.max_update_terms.min(limits.max_update_terms);
        self.system
            .validate_friction_geometry_config(FrictionGeometryConfig { max_update_terms })
            .map_err(RenderError::Friction)?;
        state.limits.max_update_terms = max_update_terms;
        state.fiber_support = Some(FiberSupportState {
            limits,
            participation,
            force_start,
            mass_kg,
            last_frame: None,
            accepted: FiberSupportKinematics {
                velocity_m_s: initial_velocity_m_s,
                displacement_m: 0.0,
                station_fraction,
            },
        });
        Ok(self)
    }
    /// Last accepted travel/speed/material crossing; refusal preserves it.
    #[must_use]
    pub fn fiber_support_kinematics(&self) -> Option<FiberSupportKinematics> {
        Some(self.bow_stations.as_ref()?.fiber_support.as_ref()?.accepted)
    }
    #[must_use]
    pub fn fiber_support_motion_limits(&self) -> Option<FiberSupportMotionLimits> {
        Some(self.bow_stations.as_ref()?.fiber_support.as_ref()?.limits)
    }
    /// Last jointly accepted support work and laboratory mechanical energy.
    #[must_use]
    pub fn fiber_support_frame(&self) -> Option<&FiberSupportFrame> {
        self.bow_stations
            .as_ref()?
            .fiber_support
            .as_ref()?
            .last_frame
            .as_ref()
    }
}
impl FiberSupportState {
    fn moments(
        &self,
        component_start: usize,
        cells: usize,
        system: &MultiFrictionModalSystem,
    ) -> (f64, f64) {
        system.components()[component_start..component_start + cells]
            .iter()
            .flat_map(|component| component.states().iter().zip(&self.participation))
            .fold((0.0, 0.0), |(q, v), (state, gamma)| {
                (
                    q + gamma * state.displacement_m_sqrt_kg,
                    v + gamma * state.velocity_m_sqrt_kg_per_s,
                )
            })
    }
    pub(crate) fn prepare(
        &self,
        source: &MultiBowStringBankSource,
        system: &MultiFrictionModalSystem,
        external: &[f64],
        drives: &[FrictionDrive],
        gate: &CancelGate,
    ) -> Result<FiberSupportStep, ModalFrictionError> {
        let fiber = admit_sources(source)?;
        let cells = fiber.model().resolution().cells;
        if drives.len() != cells || external.len() != system.mode_count() || drives.is_empty() {
            return Err(error(
                "fiber motion needs complete force and ribbon drive vectors",
            ));
        }
        let velocity = drives[0].speed_m_s;
        if !velocity.is_finite()
            || velocity.abs() > self.limits.maximum_speed_m_s
            || drives.iter().any(|d| d.speed_m_s != velocity)
        {
            return Err(error(
                "all fiber cells must share one admitted support velocity",
            ));
        }
        let dt = system.sample_period_s();
        let previous = self.accepted.velocity_m_s;
        let acceleration = (velocity - previous) / dt;
        let mean_velocity = 0.5 * previous + 0.5 * velocity;
        if !acceleration.is_finite()
            || acceleration.abs() > self.limits.maximum_acceleration_m_s2
            || dt * previous.abs().max(velocity.abs()) > self.limits.maximum_step_m
        {
            return Err(error(
                "fiber support exceeds its acceleration or travel limit",
            ));
        }
        let length = fiber.model().geometry().length_m;
        let displacement = self.accepted.displacement_m + dt * mean_velocity;
        let station = fiber.station_fraction() - displacement / length;
        // t=dt/2 under constant acceleration, not mean endpoint position.
        let midpoint_displacement =
            self.accepted.displacement_m + dt * (0.375 * previous + 0.125 * velocity);
        let midpoint_station = fiber.station_fraction() - midpoint_displacement / length;
        let interior = |s: f64| s.is_finite() && 0.0 < s && s < 1.0;
        if !interior(station) || !interior(midpoint_station) {
            return Err(error("fiber stroke crosses an original supported endpoint"));
        }
        // A reversal may leave the span and return inside the same interval.
        if previous != 0.0
            && velocity != 0.0
            && previous.is_sign_positive() != velocity.is_sign_positive()
        {
            let turning_time = -previous / acceleration;
            let turning = self.accepted.displacement_m + 0.5 * previous * turning_time;
            if !interior(fiber.station_fraction() - turning / length) {
                return Err(error("fiber reversal crosses the supported endpoint"));
            }
        }
        let shape = fiber
            .model()
            .displacement_map(midpoint_station, FiberDirection::Axial)
            .map_err(|_| error("moving fiber projection is unrepresentable"))?;
        let start_shape = fiber
            .model()
            .displacement_map(self.accepted.station_fraction, FiberDirection::Axial)
            .map_err(|_| error("initial fiber endpoint projection is unrepresentable"))?;
        let end_shape = fiber
            .model()
            .displacement_map(station, FiberDirection::Axial)
            .map_err(|_| error("final fiber endpoint projection is unrepresentable"))?;
        let mut forces = external.to_vec();
        for cell in 0..cells {
            if gate.is_requested() {
                return Err(ModalCouplingError::Cancelled.into());
            }
            let start = self.force_start + cell * self.participation.len();
            for (force, &mass) in forces[start..start + self.participation.len()]
                .iter_mut()
                .zip(&self.participation)
            {
                *force -= acceleration * mass;
                if !force.is_finite() {
                    return Err(error("fiber inertial load is unrepresentable"));
                }
            }
        }
        let drives = drives
            .iter()
            .map(|drive| FrictionDrive {
                speed_m_s: mean_velocity,
                normal_force_n: drive.normal_force_n,
            })
            .collect();
        let (previous_displacement_moment, momentum) =
            self.moments(fiber.component_start(), cells, system);
        let previous_kinetic_correction_j =
            previous * momentum + 0.5 * self.mass_kg * previous * previous;
        Ok(FiberSupportStep {
            accepted: FiberSupportKinematics {
                velocity_m_s: velocity,
                displacement_m: displacement,
                station_fraction: station,
            },
            component_start: fiber.component_start(),
            shape,
            start_shape,
            end_shape,
            start_speed_m_s: previous,
            acceleration_m_s2: acceleration,
            previous_displacement_moment,
            previous_kinetic_correction_j,
            forces,
            drives,
        })
    }
    pub(crate) fn accept(
        &mut self,
        step: FiberSupportStep,
        system: &MultiFrictionModalSystem,
        frame: &MultiFrictionFrame,
    ) {
        let (displacement, momentum) =
            self.moments(step.component_start, step.drives.len(), system);
        let velocity = step.accepted.velocity_m_s;
        let correction = velocity * momentum + 0.5 * self.mass_kg * velocity * velocity;
        let inertial_work =
            -step.acceleration_m_s2 * (displacement - step.previous_displacement_moment);
        // Material-fiber banks have no lumped series-drive storage; retained
        // fiber modes already represent that compliance. Moving shear/thermal
        // patches remain inadmissible under this owner's transport restrictions.
        let laboratory_energy = frame.network_energy_j
            + frame.obstacle_energy_j
            + frame.obstacle_residual_energy_j
            + frame.shear_layer_energy_j
            + correction;
        self.last_frame = Some(FiberSupportFrame {
            sample: frame.sample,
            velocity_m_s: velocity,
            translation_m: step.accepted.displacement_m,
            acceleration_m_s2: step.acceleration_m_s2,
            inertial_coordinate_work_j: inertial_work,
            laboratory_kinetic_correction_j: correction,
            support_work_j: frame.drive_work_j + inertial_work + correction
                - step.previous_kinetic_correction_j,
            laboratory_mechanical_energy_j: laboratory_energy,
        });
        self.accepted = step.accepted;
    }
}
