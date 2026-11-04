//! Pressure-driven motion on the actual thermal valve modes.
//!
//! The existing modal time owner advances each fixed-source interval. At a
//! source update the generic mass-inner-product transfer preserves physical
//! displacement and momentum, subject to explicit projection/work limits.
//! No mode-number matching, amplitude normalization, new integrator or thermal
//! feedback is introduced. The pressure observation is an explicitly supplied
//! narrow-band transfer; acoustic radiation reaction is not inferred from it.

use super::ThermalPressureLoadedValve;
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::time::{
    ModalAcousticFrame, ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget,
    ModalAcousticTimeError, ModalAcousticTimeModel,
};
use fs_modal::transfer::{
    ModalTransfer, ModalTransferControls, ModalTransferError, transfer_modal_state,
};

/// Original source, modal-runtime or transfer refusal. No partial state is
/// published on a failed step or source update.
#[derive(Debug)]
pub enum ThermalValveMotionError {
    Invalid(&'static str),
    Modal(ModalAcousticTimeError),
    Transfer(ModalTransferError),
    Cancelled,
}
impl core::fmt::Display for ThermalValveMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "thermal valve motion: {what}"),
            Self::Modal(error) => write!(f, "thermal valve motion: {error}"),
            Self::Transfer(error) => write!(f, "thermal valve motion: {error}"),
            Self::Cancelled => write!(f, "thermal valve motion cancelled"),
        }
    }
}
impl std::error::Error for ThermalValveMotionError {}
impl From<ModalAcousticTimeError> for ThermalValveMotionError {
    fn from(error: ModalAcousticTimeError) -> Self {
        Self::Modal(error)
    }
}
impl From<ModalTransferError> for ThermalValveMotionError {
    fn from(error: ModalTransferError) -> Self {
        Self::Transfer(error)
    }
}
fn poll(cx: &Cx<'_>) -> Result<(), ThermalValveMotionError> {
    cx.checkpoint().map_err(|_| ThermalValveMotionError::Cancelled)
}

/// One source-bound accepted pressure step on the original sample clock.
#[derive(Clone, Debug)]
pub struct ThermalValveMotionFrame {
    pub step: u64,
    pub time_s: f64,
    /// Signed pressure-port swept volume [m³], work-conjugate to the applied
    /// pressure. This excludes any change of thermal equilibrium at a rebind.
    pub swept_volume_change_m3: f64,
    /// Existing modal energy, applied pressure work, loss and observer pressure.
    pub modal: ModalAcousticFrame,
}

/// Actual reference-node displacement and velocity of the represented motion.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalValveNodalMotion {
    /// Equilibrium displacement plus the retained vibration field [m].
    pub displacement_m: Vec<[f64; 3]>,
    pub velocity_m_s: Vec<[f64; 3]>,
}

/// A thermal reduction and its continuing physical modal state.
///
/// Source refreshes are instantaneous parameter updates at the accepted sample
/// boundary. Their signed incremental modal-storage change is reported and
/// bounded, not hidden as damping or deducted from caloric heat. Time-varying
/// interpolation and reciprocal thermoelastic heat require a separate owner.
#[derive(Clone, Debug)]
pub struct ThermalValveModalMotion {
    valve: ThermalPressureLoadedValve,
    model: ModalAcousticTimeModel,
    sample_rate_hz: u32,
    accepted_steps: u64,
}

impl ThermalValveModalMotion {
    /// Start from explicit mass-normalized coordinates on this exact source
    /// basis. Supply one pressure transfer per retained mode, including explicit
    /// zero for unobserved motion; no observer or structural loss is invented.
    pub fn new(
        valve: ThermalPressureLoadedValve,
        sample_rate_hz: u32,
        initial: &[ModalAcousticState],
        pressure_transfers: &[C64],
        budget: ModalAcousticTimeBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, ThermalValveMotionError> {
        poll(cx)?;
        let mut model = modal_model(&valve, sample_rate_hz, pressure_transfers, budget, cx)?;
        model.restore_states(initial)?;
        model.observer_pressure_with_transfers(pressure_transfers)?;
        poll(cx)?;
        Ok(Self {
            valve,
            model,
            sample_rate_hz,
            accepted_steps: 0,
        })
    }

    pub const fn valve(&self) -> &ThermalPressureLoadedValve {
        &self.valve
    }
    pub const fn model(&self) -> &ModalAcousticTimeModel {
        &self.model
    }
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }
    pub fn time_s(&self) -> f64 {
        self.accepted_steps as f64 * self.model.sample_period_s()
    }

    /// Apply the original reference pressure footprint for one existing exact
    /// ZOH step. Positive pressure acts toward opening. Its transpose supplies
    /// swept volume, so modal input work equals pressure times swept volume.
    /// The pressure is prescribed here; this is not a bore-feedback solver.
    pub fn step_pressure(
        &mut self,
        pressure_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<ThermalValveMotionFrame, ThermalValveMotionError> {
        poll(cx)?;
        let next = self.accepted_steps.checked_add(1)
            .filter(|&n| n <= (1_u64 << 53))
            .ok_or(ThermalValveMotionError::Invalid("sample ordinal is not exactly representable"))?;
        let time_s = next as f64 * self.model.sample_period_s();
        if !pressure_pa.is_finite() || !time_s.is_finite() || time_s <= self.time_s() {
            return Err(ThermalValveMotionError::Invalid("finite pressure and a distinct sample clock required"));
        }
        let load = self.valve.valve().mass_normalized_pressure_load();
        let forces: Vec<_> = load.iter().map(|&area| area * pressure_pa).collect();
        let mut candidate = self.model.clone();
        let modal = candidate.step(&forces)?;
        let swept_volume_change_m3: f64 = load.iter().zip(self.model.states()).zip(candidate.states())
            .map(|((&area, old), state)| area * (state.displacement_m_sqrt_kg - old.displacement_m_sqrt_kg)).sum();
        if !swept_volume_change_m3.is_finite() {
            return Err(ThermalValveMotionError::Invalid("swept pressure volume is not finite"));
        }
        poll(cx)?;
        self.model = candidate;
        self.accepted_steps = next;
        Ok(ThermalValveMotionFrame { step: next, time_s, swept_volume_change_m3, modal })
    }

    /// Reconstruct the source mesh field, including the actual compatible
    /// equilibrium. Fixed DOFs retain their admitted zero vibration velocity.
    pub fn nodal_motion(&self, cx: &Cx<'_>)
        -> Result<ThermalValveNodalMotion, ThermalValveMotionError>
    {
        poll(cx)?;
        let mut displacement_m = self.valve.equilibrium().displacement_m().to_vec();
        let mut velocity_m_s = vec![[0.0; 3]; displacement_m.len()];
        let free = &self.valve.equilibrium().pencil().free_dofs;
        for (mode, state) in self.valve.valve().source_modes().iter().zip(self.model.states()) {
            poll(cx)?;
            for (index, (&shape, &dof)) in mode.phi.iter().zip(free).enumerate() {
                if index % 256 == 0 { poll(cx)?; }
                displacement_m[dof / 3][dof % 3] += shape * state.displacement_m_sqrt_kg;
                velocity_m_s[dof / 3][dof % 3] += shape * state.velocity_m_sqrt_kg_per_s;
            }
        }
        if !displacement_m.iter().chain(&velocity_m_s).flatten().all(|x| x.is_finite()) {
            return Err(ThermalValveMotionError::Invalid("nodal motion is not finite"));
        }
        Ok(ThermalValveNodalMotion { displacement_m, velocity_m_s })
    }

    /// Refresh the reduction while retaining physical displacement, conjugate
    /// momentum and the exact sample clock. A changed equilibrium contributes
    /// its displacement offset before projection; it never resets vibration.
    ///
    /// Reference mesh, frame, discretization, constrained DOFs and consistent
    /// mass must be identical. New source modes may reorder, rotate or change
    /// count; the generic transfer admits physical residuals and energy
    /// explicitly. A new pressure observation is supplied for that new basis.
    /// Refusal leaves this accepted source, state and clock unchanged.
    pub fn rebind_candidate(
        &self,
        next: ThermalPressureLoadedValve,
        pressure_transfers: &[C64],
        controls: ModalTransferControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, ModalTransfer), ThermalValveMotionError> {
        poll(cx)?;
        let old_card = self.valve.valve().card();
        let new_card = next.valve().card();
        let old_pencil = self.valve.equilibrium().pencil();
        let new_pencil = next.equilibrium().pencil();
        if old_card.mesh_digest != new_card.mesh_digest
            || self.valve.mesh_frame() != next.mesh_frame()
            || old_card.discretization != new_card.discretization
            || old_pencil.free_dofs != new_pencil.free_dofs
            || old_pencil.mass != new_pencil.mass
            || self.valve.equilibrium().displacement_m().len() != next.equilibrium().displacement_m().len()
        {
            return Err(ThermalValveMotionError::Invalid("source transfer requires the same reference mesh, frame, supports and consistent mass"));
        }
        let offset: Vec<_> = old_pencil.free_dofs.iter().map(|&dof| {
            self.valve.equilibrium().displacement_m()[dof / 3][dof % 3]
                - next.equilibrium().displacement_m()[dof / 3][dof % 3]
        }).collect();
        let transfer = transfer_modal_state(
            &old_pencil.mass,
            &new_pencil.stiffness,
            self.valve.valve().source_modes(),
            next.valve().source_modes(),
            self.model.states(),
            &offset,
            controls,
            cx,
        )?;
        let mut model = modal_model(&next, self.sample_rate_hz, pressure_transfers,
            self.model.budget(), cx)?;
        model.restore_states(&transfer.states)?;
        model.observer_pressure_with_transfers(pressure_transfers)?;
        poll(cx)?;
        Ok((Self { valve: next, model, sample_rate_hz: self.sample_rate_hz,
            accepted_steps: self.accepted_steps }, transfer))
    }
}

fn modal_model(
    valve: &ThermalPressureLoadedValve,
    sample_rate_hz: u32,
    transfers: &[C64],
    budget: ModalAcousticTimeBudget,
    cx: &Cx<'_>,
) -> Result<ModalAcousticTimeModel, ThermalValveMotionError> {
    let source = valve.valve();
    if transfers.len() != source.source_modes().len() {
        return Err(ThermalValveMotionError::Invalid("one pressure observation per actual retained source mode required"));
    }
    let mut modes = Vec::with_capacity(transfers.len());
    for ((pair, mode), &pressure_per_modal_velocity) in source.source_modes().iter()
        .zip(&source.card().modes).zip(transfers)
    {
        poll(cx)?;
        modes.push(ModalAcousticMode {
            angular_frequency_rad_s: pair.lambda.sqrt(),
            damping_ratio: mode.damping_n_s_m / (2.0 * (mode.mass_kg * mode.stiffness_n_m).sqrt()),
            pressure_per_modal_velocity,
        });
    }
    Ok(ModalAcousticTimeModel::try_new(sample_rate_hz, modes, budget)?)
}
