//! Transactional physical observations of the original exact modal time step.
use super::{FrameVibrationBasis, FrameVibrationError, checked_sum, poll};
use crate::modal_acoustic_time::{
    ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
    ModalAcousticWorkspace,
};

/// Endpoint physical limits, separate from the original mass-normalized limits.
/// Each is positive and finite. Translation bounds are componentwise, not norms.
#[derive(Clone, Copy, Debug)]
pub struct FrameResponseLimits {
    /// Maximum retained displacement increment in each translation coordinate [m].
    pub maximum_incremental_translation_m: f64,
    /// Maximum retained angular increment [rad].
    pub maximum_incremental_rotation_rad: f64,
    /// Maximum relative translational velocity component [m/s].
    pub maximum_translation_speed_m_s: f64,
    /// Maximum angular velocity [rad/s].
    pub maximum_angular_speed_rad_s: f64,
    /// Maximum absolute translational acceleration component [m/s^2].
    pub maximum_absolute_translation_acceleration_m_s2: f64,
    /// Maximum angular acceleration [rad/s^2].
    pub maximum_angular_acceleration_rad_s2: f64,
}

/// Only accepted endpoints are published. Arrays retain node-major [x,y,theta]
/// ordering, with translational SI units and radians for angular coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameVibrationResponse {
    /// Accepted fixed-rate steps, excluding construction.
    pub sample: u64,
    /// Elapsed physical time, derived from sample count rather than repeated addition.
    pub time_s: f64,
    /// Frozen preload displacement plus the retained modal increment.
    pub displacement: Vec<f64>,
    /// Retained modal increment, without subtracting two large preload values.
    pub incremental_displacement: Vec<f64>,
    /// Relative to the translating base; its own velocity is not integrated.
    pub relative_velocity: Vec<f64>,
    /// Equation-derived acceleration under the just-held loads, not a finite difference.
    pub relative_acceleration: Vec<f64>,
    /// Relative acceleration plus the held base acceleration on x/y entries.
    /// Fixed coordinates therefore report the actual base acceleration.
    pub absolute_acceleration: Vec<f64>,
    /// Incremental quadratic modal storage only; excludes frozen preload energy.
    pub retained_energy_j: f64,
    /// Relative-coordinate external work for THIS sample from the shared owner.
    pub input_work_j: f64,
    /// Held nodal forces/moments times retained physical displacement increments.
    pub nodal_work_j: f64,
    /// -Gamma*a_base times modal increments: not absolute foundation work.
    pub base_inertial_work_j: f64,
    /// Shared owner's work-minus-energy loss; tiny negative roundoff is retained.
    pub viscous_dissipation_j: f64,
    /// Shared modal loss roundoff allowance [J].
    pub dissipation_roundoff_tolerance_j: f64,
    /// Independent physical-work sum minus the shared owner's generalized work.
    pub work_projection_residual_j: f64,
    /// Floating-point work-conjugacy gate used on this accepted sample [J].
    pub work_projection_tolerance_j: f64,
}
impl FrameVibrationResponse {
    fn zero(n: usize) -> Self {
        Self {
            sample: 0, time_s: 0.0, displacement: vec![0.0; n],
            incremental_displacement: vec![0.0; n], relative_velocity: vec![0.0; n],
            relative_acceleration: vec![0.0; n], absolute_acceleration: vec![0.0; n],
            retained_energy_j: 0.0, input_work_j: 0.0, nodal_work_j: 0.0,
            base_inertial_work_j: 0.0, viscous_dissipation_j: 0.0,
            dissipation_roundoff_tolerance_j: 0.0, work_projection_residual_j: 0.0,
            work_projection_tolerance_j: 0.0,
        }
    }
}

/// Fixed-rate small-vibration history without a new time integrator. All step
/// buffers and both model state arrays are allocated during construction; the
/// exact full-sample transition cache is shared by the accepted and trial model.
/// Cancellation is checked around the bounded modal step and in physical loops.
/// It cannot interrupt an individual original scalar-mode propagation loop.
///
/// The original frame/material history is not modified. For nonlinear contact
/// or multibody feedback, use the same basis's `coupled_component` and attachment
/// maps with the existing `CoupledModalSystem`. That complete binding also
/// carries source section viscosity; this independent-mode wrapper refuses it
/// through `time_model` rather than silently discarding its modal cross terms.
pub struct FrameVibrationRuntime {
    basis: FrameVibrationBasis,
    accepted: ModalAcousticTimeModel,
    trial: ModalAcousticTimeModel,
    workspace: ModalAcousticWorkspace,
    forces: Vec<f64>,
    modal_acceleration: Vec<f64>,
    displacement_change: Vec<f64>,
    response: FrameVibrationResponse,
    staged: FrameVibrationResponse,
    limits: FrameResponseLimits,
    sample_rate_hz: u32,
}
impl FrameVibrationRuntime {
    /// Explicit initial coordinates in THIS retained basis. The initial physical
    /// response uses zero incremental load and stationary base. Arbitrary nodal
    /// velocities or an old nonlinear frame history are not silently imported.
    pub fn new(
        basis: FrameVibrationBasis, sample_rate_hz: u32, damping_ratios: &[f64],
        modal_budget: ModalAcousticTimeBudget, limits: FrameResponseLimits,
        initial_states: &[ModalAcousticState], mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FrameVibrationError> {
        poll(&mut cancelled)?;
        let ceilings = [limits.maximum_incremental_translation_m,
            limits.maximum_incremental_rotation_rad, limits.maximum_translation_speed_m_s,
            limits.maximum_angular_speed_rad_s,
            limits.maximum_absolute_translation_acceleration_m_s2,
            limits.maximum_angular_acceleration_rad_s2];
        if ceilings.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(FrameVibrationError::Invalid("positive finite physical response limits required"));
        }
        let mut accepted = basis.time_model(sample_rate_hz, damping_ratios, modal_budget)?;
        accepted.restore_states(initial_states)?;
        let forces = vec![0.0; basis.modes().len()];
        let mut modal_acceleration = forces.clone();
        let mut response = FrameVibrationResponse::zero(basis.ndof());
        fill_response(&basis, &accepted, &forces, [0.0; 2], limits,
            &mut modal_acceleration, &mut response, &mut cancelled)?;
        response.retained_energy_j = checked_sum(accepted.modes().iter().zip(accepted.states())
            .map(|(m,s)| 0.5*(s.velocity_m_sqrt_kg_per_s.powi(2)
                +(m.angular_frequency_rad_s*s.displacement_m_sqrt_kg).powi(2))))?;
        let trial = accepted.clone();
        let workspace = ModalAcousticWorkspace::new(&accepted);
        let staged = response.clone();
        let displacement_change = vec![0.0; basis.ndof()];
        poll(&mut cancelled)?;
        Ok(Self { basis, accepted, trial, workspace, forces, modal_acceleration,
            displacement_change, response, staged, limits, sample_rate_hz })
    }

    /// Original immutable retained structural basis.
    pub fn basis(&self) -> &FrameVibrationBasis { &self.basis }
    /// Last accepted mass-normalized state, in the retained source order.
    pub fn states(&self) -> &[ModalAcousticState] { self.accepted.states() }
    /// Last accepted physical endpoint; trial diagnostics are never exposed.
    pub fn response(&self) -> &FrameVibrationResponse { &self.response }

    /// Advance ONE fixed sample under held incremental physical loads and held
    /// two-axis uniform base acceleration. No interpolation or pulse smearing is
    /// invented. Endpoint acceleration uses the same held force and damping law.
    ///
    /// A late physical-limit, projection-work or cancellation refusal leaves
    /// accepted modal state, physical response and sample clock unchanged. Trial
    /// scratch may change. No construction, resizing or cloning occurs here.
    pub fn step(
        &mut self, nodal_force_and_moment: &[f64], base_acceleration_m_s2: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&FrameVibrationResponse, FrameVibrationError> {
        poll(&mut cancelled)?;
        let sample = self.response.sample.checked_add(1)
            .ok_or(FrameVibrationError::Invalid("frame sample counter exhausted"))?;
        let time_s = sample as f64 / f64::from(self.sample_rate_hz);
        if !time_s.is_finite() || time_s <= self.response.time_s {
            return Err(FrameVibrationError::Invalid("frame sample time is no longer representable"));
        }
        self.basis.project_force_into(nodal_force_and_moment, base_acceleration_m_s2,
            &mut self.forces, &mut cancelled)?;
        self.trial.restore_states(self.accepted.states())?;
        self.trial.step_into(&self.forces, &mut self.workspace)?;
        poll(&mut cancelled)?;
        fill_response(&self.basis, &self.trial, &self.forces, base_acceleration_m_s2,
            self.limits, &mut self.modal_acceleration, &mut self.staged, &mut cancelled)?;
        for (i, change) in self.displacement_change.iter_mut().enumerate() {
            if i % 128 == 0 { poll(&mut cancelled)?; }
            *change = checked_sum(self.basis.modes().iter()
                .zip(self.accepted.states()).zip(self.trial.states())
                .map(|((m,old),new)| m.displacement[i]
                    *(new.displacement_m_sqrt_kg-old.displacement_m_sqrt_kg)))?;
        }
        let nodal_work = checked_sum(nodal_force_and_moment.iter()
            .zip(&self.displacement_change).map(|(f,du)| f*du))?;
        let base_work_terms = || self.basis.modes().iter()
            .zip(self.accepted.states()).zip(self.trial.states())
            .map(|((m,old),new)| {
                let force = -m.base_participation[0]*base_acceleration_m_s2[0]
                    -m.base_participation[1]*base_acceleration_m_s2[1];
                force*(new.displacement_m_sqrt_kg-old.displacement_m_sqrt_kg)
            });
        let base_work = checked_sum(base_work_terms())?;
        let absolute_work = checked_sum(nodal_force_and_moment.iter()
            .zip(&self.displacement_change).map(|(f,du)| (f*du).abs())
            .chain(base_work_terms().map(f64::abs)))?;
        let modal = self.workspace.frame();
        let residual = checked_sum([nodal_work, base_work, -modal.input_work_j].into_iter())?;
        let scale = absolute_work.max(modal.input_work_j.abs())
            .max(modal.total_modal_energy_j).max(self.response.retained_energy_j)
            .max(f64::MIN_POSITIVE);
        let tolerance = 1024.0*f64::EPSILON
            *(self.basis.ndof()+self.basis.modes().len()) as f64*scale;
        if !tolerance.is_finite() || residual.abs() > tolerance {
            return Err(FrameVibrationError::Invalid("physical and modal work projections disagree"));
        }
        self.staged.sample = sample;
        self.staged.time_s = time_s;
        self.staged.retained_energy_j = modal.total_modal_energy_j;
        self.staged.input_work_j = modal.input_work_j;
        self.staged.nodal_work_j = nodal_work;
        self.staged.base_inertial_work_j = base_work;
        self.staged.viscous_dissipation_j = modal.viscous_dissipation_j;
        self.staged.dissipation_roundoff_tolerance_j = modal.dissipation_roundoff_tolerance_j;
        self.staged.work_projection_residual_j = residual;
        self.staged.work_projection_tolerance_j = tolerance;
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.accepted, &mut self.trial);
        core::mem::swap(&mut self.response, &mut self.staged);
        Ok(&self.response)
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_response(
    basis: &FrameVibrationBasis, model: &ModalAcousticTimeModel, forces: &[f64],
    base: [f64; 2], limits: FrameResponseLimits, modal_acceleration: &mut [f64],
    response: &mut FrameVibrationResponse, cancelled: &mut impl FnMut() -> bool,
) -> Result<(), FrameVibrationError> {
    for (i, ((mode, state), acceleration)) in model.modes().iter().zip(model.states())
        .zip(modal_acceleration.iter_mut()).enumerate() {
        if i % 128 == 0 { poll(cancelled)?; }
        *acceleration = checked_sum([forces[i],
            -mode.angular_frequency_rad_s*mode.angular_frequency_rad_s*state.displacement_m_sqrt_kg,
            -2.0*mode.damping_ratio*mode.angular_frequency_rad_s*state.velocity_m_sqrt_kg_per_s,
        ].into_iter())?;
    }
    basis.reconstruct_into(model.states(), &mut response.displacement,
        &mut response.relative_velocity, &mut *cancelled)?;
    for i in 0..basis.ndof() {
        if i % 128 == 0 { poll(cancelled)?; }
        let increment = checked_sum(basis.modes().iter().zip(model.states())
            .map(|(m,s)| m.displacement[i]*s.displacement_m_sqrt_kg))?;
        let acceleration = checked_sum(basis.modes().iter().zip(modal_acceleration.iter())
            .map(|(m,a)| m.displacement[i]*a))?;
        let axis = i%3;
        let absolute = checked_sum([acceleration, if axis < 2 { base[axis] } else { 0.0 }].into_iter())?;
        response.incremental_displacement[i] = increment;
        response.relative_acceleration[i] = acceleration;
        response.absolute_acceleration[i] = absolute;
        let bounds = if axis < 2 {
            [limits.maximum_incremental_translation_m, limits.maximum_translation_speed_m_s,
                limits.maximum_absolute_translation_acceleration_m_s2]
        } else {
            [limits.maximum_incremental_rotation_rad, limits.maximum_angular_speed_rad_s,
                limits.maximum_angular_acceleration_rad_s2]
        };
        for (value, limit) in [increment, response.relative_velocity[i], absolute].into_iter().zip(bounds) {
            if value.abs() > limit {
                return Err(FrameVibrationError::Invalid("physical frame displacement, velocity or acceleration limit exceeded"));
            }
        }
    }
    poll(cancelled)
}
