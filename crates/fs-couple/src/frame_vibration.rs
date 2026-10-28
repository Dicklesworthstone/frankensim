//! Physical frame loads and observations on the shared modal/time/coupling owners.
//!
//! The original `Frame2d` owns the signed tangent, consistent reference mass,
//! condensation and eigenvectors. This adapter keeps their normalization:
//! `u = u_ref + Phi*q`, `Q = Phi^T*f - Gamma*a_base`. Thus nodal forces and
//! moments do the same work as their generalized forces. There is no second
//! eigensolver, time integrator, invented radiation transfer or damping fit.
//!
//! This is small vibration about a frozen, externally maintained reference
//! state, not nonlinear continuation of its material/geometry history. Loads
//! are increments about that reference; do not apply its gravity preload twice.
//! Omitted modes and supported inertia remain explicit in the mass inventory.
//! Outputs lie in the retained modal subspace: no residual static flexibility
//! or instantaneous load-dependent massless-coordinate feedthrough is inferred.

use fs_math::c64::C64;
use fs_solid::frame2d::{Frame2d, FrameMode, FrameModeStatus};
use fs_solid::SolidError;

use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget,
    ModalAcousticTimeError, ModalAcousticTimeModel,
};
use crate::render::schedule::force::coupled::ModalAttachment;

/// Admission before the original dense frame assembly/eigensolve allocates.
#[derive(Clone, Copy, Debug)]
pub struct FrameVibrationBudget {
    /// Full physical coordinates, including fixed and massless coordinates.
    pub max_frame_dofs: usize,
    /// Full `ndof*ndof` entries, not merely the retained modal output.
    /// The native frame/eigensolver has its own additional admission limits.
    pub max_dense_entries: usize,
}

/// A refused frame-to-modal binding or physical projection.
#[derive(Debug)]
pub enum FrameVibrationError {
    /// Invalid dimensions, physical domain, budget or unrepresentable result.
    Invalid(&'static str),
    /// Original structural assembly or eigensolver refusal.
    Solid(SolidError),
    /// Original shared time-domain owner refusal.
    Modal(ModalAcousticTimeError),
    /// No unstable mode may be hidden by retaining a smaller positive subset.
    UnstableMode { mode: usize, eigenvalue_per_s2: f64 },
    /// A near-zero mode cannot be silently converted into a positive oscillator.
    UnresolvedMode { mode: usize, eigenvalue_per_s2: f64 },
    /// Cancellation before publication; output scratch may have been overwritten.
    Cancelled,
}
impl core::fmt::Display for FrameVibrationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "frame vibration: {self:?}")
    }
}
impl std::error::Error for FrameVibrationError {}
impl From<SolidError> for FrameVibrationError {
    fn from(error: SolidError) -> Self { Self::Solid(error) }
}
impl From<ModalAcousticTimeError> for FrameVibrationError {
    fn from(error: ModalAcousticTimeError) -> Self { Self::Modal(error) }
}

/// An existing modal mechanical port plus its frozen physical reference datum.
#[derive(Clone, Debug)]
pub struct FrameTranslationAttachment {
    /// The original reciprocal force/displacement map, in incremental coordinates.
    pub port: ModalAttachment,
    /// Add to the port displacement to recover total reference-relative position.
    /// Account for this datum when specifying spring rest extensions/contact gaps.
    pub reference_displacement_m: f64,
}

/// Immutable source-derived basis, shared by forcing, observations and contacts.
#[derive(Clone, Debug)]
pub struct FrameVibrationBasis {
    modes: Vec<FrameMode>,
    reference_displacement: Vec<f64>,
    available_modes: usize,
    influence_mass_kg: [f64; 2],
    effective_mass_kg: [f64; 2],
}

impl FrameVibrationBasis {
    /// Extract the complete signed spectrum, refuse any unstable/unresolved
    /// mode, then explicitly retain the lowest `retained_modes` oscillators.
    /// No Nyquist-based or stability-based silent truncation occurs.
    ///
    /// Cancellation is checked around the bounded native eigensolve and per
    /// mode; this adapter cannot preempt that synchronous eigensolve itself.
    /// Existing dynamic velocity is NOT imported: the returned time model is
    /// initially at rest about `frame.displacements()`.
    pub fn from_frame(
        frame: &Frame2d,
        retained_modes: usize,
        budget: FrameVibrationBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FrameVibrationError> {
        poll(&mut cancelled)?;
        let n = frame.ndof();
        if n == 0 || n > budget.max_frame_dofs || retained_modes == 0
            || n.checked_mul(n).is_none_or(|entries| entries > budget.max_dense_entries)
        {
            return Err(FrameVibrationError::Invalid("full frame/dense preparation budget or retained mode count"));
        }
        let mut analysis = frame.modal_analysis()?;
        poll(&mut cancelled)?;
        if retained_modes > analysis.modes.len() {
            return Err(FrameVibrationError::Invalid("retained count exceeds the complete frame spectrum"));
        }
        for (index, mode) in analysis.modes.iter().enumerate() {
            poll(&mut cancelled)?;
            match mode.status() {
                FrameModeStatus::Oscillatory => {},
                FrameModeStatus::Unstable => return Err(FrameVibrationError::UnstableMode {
                    mode: index, eigenvalue_per_s2: mode.eigenvalue_per_s2,
                }),
                FrameModeStatus::Unresolved => return Err(FrameVibrationError::UnresolvedMode {
                    mode: index, eigenvalue_per_s2: mode.eigenvalue_per_s2,
                }),
            }
        }
        let available_modes = analysis.modes.len();
        analysis.modes.truncate(retained_modes);
        let effective_mass_kg = [
            checked_sum(analysis.modes.iter().map(|m| m.effective_mass_kg[0]))?,
            checked_sum(analysis.modes.iter().map(|m| m.effective_mass_kg[1]))?,
        ];
        let reference_displacement = frame.displacements().to_vec();
        if reference_displacement.iter().any(|x| !x.is_finite()) {
            return Err(FrameVibrationError::Invalid("nonfinite frozen frame displacement"));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            modes: analysis.modes, reference_displacement, available_modes,
            influence_mass_kg: analysis.base_influence_mass_kg, effective_mass_kg,
        })
    }

    /// Retained source eigenpairs, with their residual and resolution diagnostics.
    pub fn modes(&self) -> &[FrameMode] { &self.modes }
    /// Full physical coordinate count, node-major `[ux, uy, theta]`.
    pub fn ndof(&self) -> usize { self.reference_displacement.len() }
    /// Complete stable source spectrum count before explicit truncation.
    pub fn available_mode_count(&self) -> usize { self.available_modes }
    /// Frozen accepted displacement; dynamic outputs are increments about it.
    pub fn reference_displacement(&self) -> &[f64] { &self.reference_displacement }
    /// Full horizontal/vertical reference influence inventory, including support mass.
    pub fn base_influence_mass_kg(&self) -> [f64; 2] { self.influence_mass_kg }
    /// Sum of retained squared participation factors, not the complete inventory.
    pub fn retained_effective_mass_kg(&self) -> [f64; 2] { self.effective_mass_kg }

    /// Construct the ORIGINAL exact held-force modal owner, with explicitly
    /// supplied per-mode viscous damping ratios. No acoustic pressure transfer
    /// is guessed from a planar structural model. This component can also be
    /// passed directly to `CoupledModalSystem` and its existing contact wrappers.
    pub fn time_model(
        &self,
        sample_rate_hz: u32,
        damping_ratios: &[f64],
        budget: ModalAcousticTimeBudget,
    ) -> Result<ModalAcousticTimeModel, FrameVibrationError> {
        if damping_ratios.len() != self.modes.len() {
            return Err(FrameVibrationError::Invalid("one explicit damping ratio per retained mode is required"));
        }
        let modes = self.modes.iter().zip(damping_ratios).map(|(mode, &damping_ratio)| {
            ModalAcousticMode {
                angular_frequency_rad_s: mode.eigenvalue_per_s2.sqrt(),
                damping_ratio,
                pressure_per_modal_velocity: C64::ZERO,
            }
        }).collect();
        Ok(ModalAcousticTimeModel::try_new(sample_rate_hz, modes, budget)?)
    }

    /// Project incremental nodal forces [N] and moments [N m], together with
    /// simultaneous horizontal/vertical uniform base acceleration [m/s^2].
    /// `Q = Phi^T*f - Gamma_x*a_x - Gamma_y*a_y` in N/sqrt(kg).
    /// Relative-frame inertial work is not absolute actuator work.
    ///
    /// No allocation; `generalized` is scratch and may be partly overwritten
    /// on overflow/cancellation. No accepted physical state is changed.
    pub fn project_force_into(
        &self,
        nodal_force_and_moment: &[f64],
        base_acceleration_m_s2: [f64; 2],
        generalized: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), FrameVibrationError> {
        poll(&mut cancelled)?;
        if nodal_force_and_moment.len() != self.ndof() || generalized.len() != self.modes.len()
            || nodal_force_and_moment.iter().chain(base_acceleration_m_s2.iter()).any(|x| !x.is_finite())
        {
            return Err(FrameVibrationError::Invalid("finite physical loads, base acceleration and matching modal output required"));
        }
        for (mode, force) in self.modes.iter().zip(generalized) {
            poll(&mut cancelled)?;
            let applied = checked_sum(mode.displacement.iter().zip(nodal_force_and_moment).map(|(p, f)| p*f))?;
            *force = checked_sum([
                applied,
                -mode.base_participation[0]*base_acceleration_m_s2[0],
                -mode.base_participation[1]*base_acceleration_m_s2[1],
            ].into_iter())?;
        }
        poll(&mut cancelled)
    }

    /// Reconstruct total reference-relative displacement and relative velocity.
    /// Translational entries use m and m/s; angular entries use rad and rad/s.
    /// The base's own translation/velocity is deliberately not added.
    /// Output buffers are scratch on a refusal. No allocation occurs.
    pub fn reconstruct_into(
        &self,
        states: &[ModalAcousticState],
        displacement: &mut [f64],
        velocity: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), FrameVibrationError> {
        poll(&mut cancelled)?;
        if states.len() != self.modes.len() || displacement.len() != self.ndof() || velocity.len() != self.ndof()
            || states.iter().any(|s| !s.displacement_m_sqrt_kg.is_finite() || !s.velocity_m_sqrt_kg_per_s.is_finite())
        {
            return Err(FrameVibrationError::Invalid("finite modal state and complete physical output buffers required"));
        }
        for i in 0..self.ndof() {
            if i % 128 == 0 { poll(&mut cancelled)?; }
            displacement[i] = checked_sum(core::iter::once(self.reference_displacement[i]).chain(
                self.modes.iter().zip(states).map(|(mode, state)| mode.displacement[i]*state.displacement_m_sqrt_kg)))?;
            velocity[i] = checked_sum(self.modes.iter().zip(states)
                .map(|(mode, state)| mode.displacement[i]*state.velocity_m_sqrt_kg_per_s))?;
        }
        poll(&mut cancelled)
    }

    /// Mechanical attachment at an original node along an explicitly unit
    /// translation direction. The SAME shape maps contact force into modes and
    /// modal motion back to the attachment, preserving virtual work.
    /// These are incremental ports; the separate reference datum is not lost.
    pub fn translation_attachment(
        &self,
        component: usize,
        node: usize,
        direction: [f64; 2],
    ) -> Result<FrameTranslationAttachment, FrameVibrationError> {
        if node >= self.ndof()/3 || direction.iter().any(|x| !x.is_finite())
            || (direction[0].hypot(direction[1]) - 1.0).abs() > 1e-12
        {
            return Err(FrameVibrationError::Invalid("existing frame node and finite unit attachment direction required"));
        }
        let coordinate = |shape: &[f64]| checked_sum([
            direction[0]*shape[3*node], direction[1]*shape[3*node+1],
        ].into_iter());
        let shapes = self.modes.iter().map(|mode| coordinate(&mode.displacement))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(FrameTranslationAttachment {
            port: ModalAttachment { component, shapes },
            reference_displacement_m: coordinate(&self.reference_displacement)?,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), FrameVibrationError> {
    if cancelled() { Err(FrameVibrationError::Cancelled) } else { Ok(()) }
}

/// Compensated deterministic reduction for physical power-conjugate maps.
fn checked_sum(values: impl Iterator<Item = f64>) -> Result<f64, FrameVibrationError> {
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        if !value.is_finite() { return Err(FrameVibrationError::Invalid("nonfinite frame projection term")); }
        let next = sum + value;
        correction += if sum.abs() >= value.abs() { (sum-next)+value } else { (value-next)+sum };
        sum = next;
    }
    let result = sum + correction;
    if result.is_finite() { Ok(result) } else { Err(FrameVibrationError::Invalid("frame projection overflow")) }
}
