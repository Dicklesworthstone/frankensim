//! A prescribed exterior gas at the existing terminal pressure/volume-flow port.
//!
//! The bore keeps its original material, travel time and viscothermal laws.
//! Every radiating mouth and exterior observer consume this medium. This is
//! a stationary, pressure-matched acoustic interface, not an advective mixing
//! layer, thermal equilibration or live replacement of an acoustic state.

use super::{DuctError, DuctResponse, GasState, LossModel, OpenHoleFlow, Region,
    RegionalDuct, Termination, WallModel, invalid};

impl RegionalDuct {
    /// Bind the homogeneous exterior independently of the internal terminal gas.
    ///
    /// This cold construction changes no geometry or stored acoustic state.
    /// Omitting it uses the actual outlet gas at every exposed mouth and
    /// observer, including upstream holes and continuous fields. In particular,
    /// it never assigns a different exterior to each hole's local bore gas.
    /// An explicitly supplied identical gas also takes the original uniform path.
    /// The exterior has the common static pressure but may differ in temperature
    /// or composition. Closed/ideal-open end conditions do not acquire radiation.
    ///
    /// Open/vented holes propagate through their original bore-gas chimney,
    /// terminated by this exterior's radiation impedance. Pressure and volume
    /// flow meet at the original mouth, with no second end correction. Closed
    /// chimney compliance is unchanged. Existing mutual-hole restrictions remain.
    ///
    /// # Errors
    /// Incoherent gas or unequal static pressure.
    pub fn with_exterior_gas(mut self, gas: GasState) -> Result<Self, DuctError> {
        gas.validate_coherence().map_err(|_| invalid("exterior gas is not a coherent material state"))?;
        if gas.pressure != self.inlet_gas().pressure {
            return Err(invalid("exterior and duct need a common static pressure; mean flow is not modeled"));
        }
        self.exterior_gas = Some(gas);
        Ok(self)
    }

    /// Homogeneous exterior used by radiation loading and physical observation.
    /// Without an explicit binding this is the original terminal segment's gas,
    /// not the ambient fallback or the inlet medium.
    #[must_use]
    pub fn exterior_gas(&self) -> &GasState {
        self.exterior_gas.as_ref().unwrap_or_else(|| self.outlet_gas())
    }

    pub(super) fn terminal_response(&self, terminal: &Region, omega: f64,
        loss: LossModel, termination: Termination, wall: WallModel<'_>)
        -> Result<DuctResponse, DuctError> {
        fs_duct::exterior::input_impedance(&terminal.duct, &terminal.gas,
            self.exterior_gas(), omega, loss, termination, wall)
    }

    pub(super) fn terminal_response_with_holes(
        &self,
        terminal: &Region,
        omega: f64,
        loss: LossModel,
        termination: Termination,
        wall: WallModel<'_>,
    ) -> Result<(DuctResponse, Vec<OpenHoleFlow>), DuctError> {
        fs_duct::exterior::input_impedance_and_holes(&terminal.duct, &terminal.gas,
            self.exterior_gas(), omega, loss, termination, wall)
    }
}

#[cfg(test)]
mod tests;
