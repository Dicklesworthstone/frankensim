//! Original finite-gas heat/work feedback for retained material-valve mechanics.
use super::{LoadedDgPerformance, invalid};
use crate::bernoulli_aperture::dg::finite_gas::FiniteGasApertureDgStep;
use crate::bernoulli_aperture::performance::finite_gas::{
    FiniteGasHeatProgram, FiniteGasPlayback, GasHeatControl,
};
use crate::render::schedule::GestureCompileError;
use fs_exec::CancelGate;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use std::sync::Arc;

impl LoadedDgPerformance {
    /// Attach heat to the original uniform frozen-gas inventory before playback.
    ///
    /// The existing caloric and acoustic work transaction refreshes this SAME
    /// loaded valve and DG state. Every material mode, blade-force control,
    /// spatial traction and material history remains in its original owner.
    /// Powers use the mechanical clock and publish only with accepted samples.
    ///
    /// This is the closed-inventory thermal experiment with its original internal
    /// pressure sensor. Computed exterior radiation, independent branch gases,
    /// upstream tracts and fixed supply passages require explicit mass/enthalpy
    /// exchange or additional inventories and remain refused. Numerical losses
    /// are not silently deposited as heat.
    pub fn with_finite_gas_heat(
        mut self, gas: ClosedGasParcel, program: FiniteGasHeatProgram, gate: Arc<CancelGate>,
    ) -> Result<Self, GestureCompileError> {
        if self.poisoned || self.system.aperture().accepted_steps() != 0
            || self.finite_gas.is_some()
        {
            return Err(invalid("finite gas heat must bind once before loaded DG playback"));
        }
        if self.system.upstream_load().is_some() || self.system.supply_channel().is_some() {
            return Err(invalid("finite gas heat requires independent inventories for an upstream tract or fixed supply passage"));
        }
        self.finite_gas = Some(FiniteGasPlayback::new(
            gas, program, gate, self.system.load(), self.system.aperture().time_step_s(),
            self.config.samples, self.remaining_control_budgets()?,
        )?);
        Ok(self)
    }

    /// Accepted finite material inventory, including its original source model.
    #[must_use]
    pub fn finite_gas(&self) -> Option<&ClosedGasParcel> {
        self.finite_gas.as_ref().map(|state| &state.gas)
    }
    /// Authored powers and original caloric/acoustic work solve admissions.
    #[must_use]
    pub fn finite_gas_heat_program(&self) -> Option<&FiniteGasHeatProgram> {
        self.finite_gas.as_ref().map(|state| &state.program)
    }
    /// Shared cancellation gate retained throughout physical playback.
    #[must_use]
    pub fn finite_gas_cancel_gate(&self) -> Option<&Arc<CancelGate>> {
        self.finite_gas.as_ref().map(|state| &state.gate)
    }
    /// Heat controls published with accepted mechanical/flow samples.
    #[must_use]
    pub fn applied_heat_controls(&self) -> &[GasHeatControl] {
        self.finite_gas.as_ref().map_or(&[], |state| &state.program.controls[..state.next])
    }
    /// Original complete material/mechanical/acoustic/thermal work report.
    #[must_use]
    pub fn last_finite_gas_step(&self) -> Option<&FiniteGasApertureDgStep<GasCompositionStep>> {
        self.finite_gas.as_ref().and_then(|state| state.last.as_ref())
    }
}
