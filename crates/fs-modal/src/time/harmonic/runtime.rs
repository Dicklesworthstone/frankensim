//! Harmonic loading through the shared transactional mechanical runtime.

use super::super::{
    ExactZohKernel, ModalAcousticFrame, ModalAcousticTimeError, ModalAcousticTimeModel,
    ModalAcousticWorkspace, ModalStepData, stage_modal_frame,
};
use super::{PreparedHarmonicForce, invalid};
use std::sync::Arc;

/// Prepared continuous harmonic load for every retained mode of one shared
/// mechanical model. Each mode receives its own cosine/sine generalized force
/// [N/sqrt(kg)] at the same nonnegative forcing frequency [rad/s]. Coefficients
/// are interval-local, with the same `exp(-i Omega t)` convention as
/// [`super::HarmonicForceTransition`]. Preparation allocates; stepping does not.
///
/// Mode count, frequencies and damping are bound to the source model's shared
/// immutable kernel. Its clones/checkpoint restores may reuse preparation;
/// an unrelated model cannot silently reuse different physical coefficients.
/// Eight-point Gaussian power integration is confined to a resolved short
/// interval and composed to the requested duration. The positive heat factor
/// and signed work form are independent of endpoint energy subtraction.
/// Their arithmetic allowances are numerical estimates, not interval bounds.
///
/// This is one continuous forcing frequency, including DC. Cosine and sine
/// quadrature interference is retained. Heat from separately prepared tones
/// must not be added to model simultaneous frequencies: their velocities have
/// additional cross terms. No modes are discarded to meet a clock or budget.
///
/// The source model must have zero legacy narrow-band pressure transfers.
/// Use [`ModalAcousticTimeModel::without_pressure_transfers`] and a separate
/// physical broadband observer when sound is required: a transfer only at a
/// mode's natural frequency does not describe a different driven frequency
/// together with its startup transient. Returned pressure is consequently zero.
#[derive(Clone, Debug)]
pub struct PreparedModalHarmonic {
    owner: Arc<ExactZohKernel>,
    modes: Box<[PreparedHarmonicForce]>,
    duration_s: f64,
    forcing_angular_frequency_rad_s: f64,
}

impl PreparedModalHarmonic {
    /// Prepare one physical forcing frequency and positive interval duration.
    /// The model's nominal Nyquist guard admits the prescribed force as well
    /// as its already-admitted natural frequencies. A longer integration window
    /// can contain several cycles; it does not resolve observations inside it.
    ///
    /// # Errors
    /// Invalid/unrepresentable frequency or duration, forcing outside the
    /// model's admitted band, nonzero narrow-band pressure transfers, or an
    /// unrepresentable motion/work/heat coefficient.
    pub fn try_new(
        model: &ModalAcousticTimeModel,
        forcing_angular_frequency_rad_s: f64,
        duration_s: f64,
    ) -> Result<Self, ModalAcousticTimeError> {
        require_mechanical_observation(model)?;
        let maximum =
            core::f64::consts::PI * model.budget.nyquist_guard_fraction / model.sample_period_s;
        if !forcing_angular_frequency_rad_s.is_finite()
            || forcing_angular_frequency_rad_s < 0.0
            || forcing_angular_frequency_rad_s > maximum
        {
            return Err(invalid(
                "harmonic forcing must lie inside the model's nominal Nyquist guard",
            ));
        }
        let mut modes = Vec::with_capacity(model.modes.len());
        for mode in &model.modes {
            modes.push(PreparedHarmonicForce::try_new(
                mode.angular_frequency_rad_s,
                mode.damping_ratio,
                forcing_angular_frequency_rad_s,
                duration_s,
            )?);
        }
        Ok(Self {
            owner: Arc::clone(&model.transition_kernel),
            modes: modes.into_boxed_slice(),
            duration_s,
            forcing_angular_frequency_rad_s,
        })
    }

    /// Physical duration of every prepared step [s].
    #[must_use]
    pub fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Shared continuous forcing frequency [rad/s], including exact DC.
    #[must_use]
    pub fn forcing_angular_frequency_rad_s(&self) -> f64 {
        self.forcing_angular_frequency_rad_s
    }

    /// Propose one mode's endpoint for a coupled continuous observer.
    /// This pure operation uses exactly the motion used by the physical step;
    /// it does not perform the model's state/energy admission or commit state.
    ///
    /// # Errors
    /// Invalid mode index, nonfinite input, or unrepresentable motion.
    pub fn advance_mode(
        &self,
        mode_index: usize,
        state: super::super::ModalAcousticState,
        cosine_force: f64,
        sine_force: f64,
    ) -> Result<super::super::ModalAcousticState, ModalAcousticTimeError> {
        self.modes
            .get(mode_index)
            .ok_or_else(|| invalid("harmonic mode index must be admitted"))?
            .motion()
            .advance(state, cosine_force, sine_force)
    }

    /// Rotate a mode's force coefficients to the next prepared interval.
    /// Advance the caller's force phase only after a successful mechanical
    /// step. A checkpoint retains both mechanical state and force phase.
    ///
    /// # Errors
    /// Nonfinite or unrepresentable rotated force coefficients.
    pub fn advance_force_phase(
        &self,
        cosine_force: f64,
        sine_force: f64,
    ) -> Result<(f64, f64), ModalAcousticTimeError> {
        self.modes[0]
            .motion()
            .advance_force_phase(cosine_force, sine_force)
    }
}

impl ModalAcousticTimeModel {
    /// Advance every mode under continuous harmonic force, retaining physical
    /// input work and independent positive viscous heat in the ordinary frame.
    /// The existing state, energy and work/heat balance gates apply unchanged.
    /// A failure leaves every accepted modal state unchanged. Workspace scratch
    /// can change; no frame is returned for a rejected step.
    ///
    /// The force arrays are cosine/sine coefficients at the interval start.
    /// This method does not mutate external phase: after success, the caller
    /// advances it with [`PreparedModalHarmonic::advance_force_phase`].
    ///
    /// # Errors
    /// Preparation from a different model, wrong/nonfinite force arrays,
    /// nonzero narrow-band pressure transfer, incompatible workspace, physical
    /// budget crossing, or disagreement of independent work/heat and energy.
    pub fn step_harmonic_into<'w>(
        &mut self,
        prepared: &PreparedModalHarmonic,
        cosine_forces: &[f64],
        sine_forces: &[f64],
        workspace: &'w mut ModalAcousticWorkspace,
    ) -> Result<&'w ModalAcousticFrame, ModalAcousticTimeError> {
        if !Arc::ptr_eq(&prepared.owner, &self.transition_kernel) {
            return Err(invalid(
                "harmonic preparation must belong to this model or its clone",
            ));
        }
        require_mechanical_observation(self)?;
        for forces in [cosine_forces, sine_forces] {
            if forces.len() != self.modes.len() {
                return Err(ModalAcousticTimeError::ForceCountMismatch {
                    expected: self.modes.len(),
                    found: forces.len(),
                });
            }
            if forces.iter().any(|x| !x.is_finite()) {
                return Err(invalid("harmonic generalized forces must be finite"));
            }
        }
        stage_modal_frame(
            &self.modes,
            &self.states,
            self.budget,
            workspace,
            "independent harmonic work, viscous heat and modal energy disagree",
            "aggregate harmonic work, viscous heat and modal energy disagree",
            |index, _mode, state| {
                let mode = &prepared.modes[index];
                let fc = cosine_forces[index];
                let fs = sine_forces[index];
                let frame = mode.advance_with_energy(state, fc, fs)?;
                Ok(ModalStepData {
                    next: frame.state,
                    work: frame.input_work_j,
                    work_allowance: frame.input_work_roundoff_tolerance_j,
                    loss: frame.viscous_dissipation_j,
                    loss_allowance: frame.dissipation_roundoff_tolerance_j,
                })
            },
        )?;
        self.states.copy_from_slice(&workspace.candidate);
        Ok(&workspace.frame)
    }
}

pub(super) fn require_mechanical_observation(
    model: &ModalAcousticTimeModel,
) -> Result<(), ModalAcousticTimeError> {
    if model.modes.iter().any(|mode| {
        mode.pressure_per_modal_velocity.re != 0.0 || mode.pressure_per_modal_velocity.im != 0.0
    }) {
        return Err(invalid(
            "harmonic motion requires zero narrow-band pressure transfers; use a separate broadband observer",
        ));
    }
    Ok(())
}
