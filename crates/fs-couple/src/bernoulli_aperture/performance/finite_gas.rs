//! Heat controls on the existing moving-aperture/finite-gas transaction.
//!
//! Heat is physical power held on the original mechanical clock. The original
//! caloric owner pays acoustic parameter work and updates the SAME source gas,
//! wave state, material wall and terminal histories. This adds no mean mass flow,
//! condensation, wall heating, or new acoustic/thermal integration law.

use super::{AperturePerformance, CoupledAperture, GestureCompileError, invalid, receiver_error};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dg::{ApertureDg, DgApertureFrame};
use crate::bernoulli_aperture::dg::finite_gas::FiniteGasApertureDgStep;
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::dg_port::finite_gas::{FiniteGasDgError, FiniteGasDgSolve, admit_initial};
use crate::render::RenderError;
use fs_exec::CancelGate;
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use std::sync::Arc;

/// Signed external heat power, applied before the named mechanical sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasHeatControl {
    /// Source sample in the performance's half-open window.
    pub sample: u64,
    /// Watts into the finite gas; negative values prescribe cooling.
    pub power_w: f64,
}

/// An explicit finite heat history and the original two nonlinear solve budgets.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasHeatProgram {
    /// Strictly increasing samples, starting at zero. The last power is held.
    pub controls: Vec<GasHeatControl>,
    /// Original source-supported NASA-9 energy inversion.
    pub caloric: GasEnergySolve,
    /// Original reciprocal gas/acoustic work and partition-step admissions.
    pub coupling: FiniteGasDgSolve,
}

pub(super) struct FiniteGasPlayback {
    gas: ClosedGasParcel,
    pub(super) program: FiniteGasHeatProgram,
    gate: Arc<CancelGate>,
    next: usize,
    held_power_w: f64,
    last: Option<FiniteGasApertureDgStep<GasCompositionStep>>,
}

fn runtime(error: FiniteGasDgError<GasEnergyError>) -> RenderError {
    match error {
        FiniteGasDgError::Cancelled | FiniteGasDgError::Background(GasEnergyError::Cancelled) =>
            RenderError::Voice(AcousticRealizeError::Cancelled),
        FiniteGasDgError::Acoustic(error) => RenderError::Voice(error),
        error => receiver_error(error.to_string()),
    }
}

impl FiniteGasPlayback {
    fn candidate(&self, sample: u64) -> f64 {
        self.program.controls.get(self.next).filter(|event| event.sample == sample)
            .map_or(self.held_power_w, |event| event.power_w)
    }

    pub(super) fn step(&mut self, system: &mut ApertureDg, drive: TubeDrive, force_n: f64)
        -> Result<DgApertureFrame, RenderError>
    {
        let heat_j = self.candidate(system.aperture().accepted_steps())
            * system.aperture().spec().time_step_s;
        let report = system.step_with_finite_gas_heat_controlled(
            &mut self.gas, drive, force_n, heat_j, self.program.caloric,
            self.program.coupling, &mut || self.gate.is_requested(),
        ).map_err(runtime)?;
        let frame = report.coupled.clone();
        self.last = Some(report);
        Ok(frame)
    }

    pub(super) fn accept(&mut self, sample: u64) {
        if let Some(event) = self.program.controls.get(self.next).filter(|event| event.sample == sample) {
            self.held_power_w = event.power_w;
            self.next += 1;
        }
    }
}

impl AperturePerformance {
    /// Attach finite thermochemical heating to this unadvanced DG source.
    ///
    /// The inventory must match the original uniform circular main bore and its
    /// finite terminal passages. An independent upstream tract needs its own
    /// inventory and is refused here. The original source/wall/diffusion/terminal
    /// owners retain their supported domains. Computed exterior radiation and
    /// independent branch gases require mean-pressure/mass exchange and refuse
    /// this closed-inventory route. Internal pressure remains available through
    /// the original observer. Numerical and exported losses are
    /// not silently deposited as heat.
    ///
    /// The gate is retained throughout playback. A refused physical sample
    /// preserves its gas, source and acoustic state; PressureRenderer then
    /// poisons the callback in the same way as the existing fixed-gas path.
    pub fn with_finite_gas_heat(mut self, gas: ClosedGasParcel,
        program: FiniteGasHeatProgram, gate: Arc<CancelGate>)
        -> Result<Self, GestureCompileError>
    {
        if self.poisoned || self.system.aperture().accepted_steps() != 0
            || self.finite_gas.is_some() || self.finite_transport.is_some()
            || self.finite_wall.is_some()
        {
            return Err(invalid("finite gas heat must bind once before playback; combined finite gas/body heat needs a joint transaction"));
        }
        let (assignments, work) = self.remaining_control_budgets()?;
        if program.controls.is_empty() || program.controls.len() > assignments
            || program.controls.len() as u128 > u128::from(work)
            || program.controls[0].sample != 0
        {
            return Err(invalid("finite gas heat requires an initial assignment and shared control/work capacity"));
        }
        let dt = self.system.aperture().spec().time_step_s;
        let mut previous = None;
        for event in &program.controls {
            if gate.is_requested() {
                return Err(GestureCompileError::Render(RenderError::Voice(AcousticRealizeError::Cancelled)));
            }
            if event.sample >= self.config.samples || previous.is_some_and(|sample| event.sample <= sample)
                || !event.power_w.is_finite() || !(event.power_w * dt).is_finite()
                || (event.power_w != 0.0 && event.power_w * dt == 0.0)
            {
                return Err(invalid("heat assignments require increasing in-window samples and representable signed power"));
            }
            previous = Some(event.sample);
        }
        program.caloric.validate().map_err(|error|
            GestureCompileError::Render(receiver_error(error.to_string())))?;
        let temperature = gas.state().state().temperature;
        if temperature < program.caloric.minimum_temperature_k
            || temperature > program.caloric.maximum_temperature_k
        {
            return Err(invalid("finite gas caloric bracket must contain the original temperature"));
        }
        let CoupledAperture::Dg(system) = &self.system else {
            return Err(invalid("finite thermochemical heat requires the original DG aperture owner"));
        };
        if system.upstream_load().is_some() {
            return Err(invalid("one heated parcel cannot replace independent upstream and downstream inventories"));
        }
        if system.load().radiation_record().is_some() || !system.load().branches().is_empty() {
            return Err(invalid("closed-inventory heating cannot impose a mean-pressure jump on a computed exterior or independent branch; transport is required"));
        }
        if system.load().duct().gas_field().is_none_or(|field| field.phase() != ThermochemicalPhase::FrozenGas) {
            return Err(invalid("finite gas heat requires the original explicit frozen-gas field"));
        }
        admit_initial::<GasEnergyError>(system.load(), &gas, program.coupling,
            &mut || gate.is_requested()).map_err(|error| GestureCompileError::Render(runtime(error)))?;
        // Admit the selected existing loss/transfer owner without advancing or
        // publishing any history. An incompatible refresh mode fails during
        // attachment, before an output consumer can start a physical sample.
        let mut admitted = system.load().clone();
        admitted.rebind_gas_field_controlled(
            system.load().duct().gas_field().expect("admitted original field"),
            program.coupling.loss.into(),
            program.coupling.work_bracket_j[0].abs().max(program.coupling.work_bracket_j[1].abs()),
            &mut || gate.is_requested(),
        ).map_err(|error| GestureCompileError::Render(RenderError::Voice(error)))?;
        self.finite_gas = Some(FiniteGasPlayback {
            gas, program, gate, next: 0, held_power_w: 0.0, last: None,
        });
        Ok(self)
    }

    /// The actual accepted inventory, including its conserved mass and source model.
    #[must_use]
    pub fn finite_gas(&self) -> Option<&ClosedGasParcel> {
        self.finite_gas.as_ref().map(|state| &state.gas).or_else(|| {
            let id = self.finite_gas_transport_program()?.source_id;
            self.finite_gas_parcels()?.get(&id)
        })
    }
    /// Original authored powers and solver admissions.
    #[must_use]
    pub fn finite_gas_heat_program(&self) -> Option<&FiniteGasHeatProgram> {
        self.finite_gas.as_ref().map(|state| &state.program)
            .or_else(|| self.finite_gas_transport_program().map(|program| &program.source))
    }
    /// Shared cancellation gate retained by the finite-gas runtime.
    #[must_use]
    pub fn finite_gas_cancel_gate(&self) -> Option<&Arc<CancelGate>> {
        self.finite_gas.as_ref().map(|state| &state.gate)
            .or_else(|| self.finite_gas_transport_cancel_gate())
    }
    /// Heat assignments published with successful pressure samples.
    #[must_use]
    pub fn applied_heat_controls(&self) -> &[GasHeatControl] {
        self.finite_gas.as_ref().map_or_else(|| self.applied_gas_transport_heat_controls(),
            |state| &state.program.controls[..state.next])
    }
    /// The original full heat-only gas/acoustic/mechanical work report.
    /// Reservoir exchange retains this source report inside
    /// [Self::last_finite_gas_transport_step] under the authored source ID.
    #[must_use]
    pub fn last_finite_gas_step(&self) -> Option<&FiniteGasApertureDgStep<GasCompositionStep>> {
        self.finite_gas.as_ref().and_then(|state| state.last.as_ref())
    }
}
