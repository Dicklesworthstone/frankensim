//! A prescribed exterior gas at the existing terminal pressure/volume-flow port.
//!
//! The bore keeps its original material, travel time and viscothermal laws.
//! Only the radiating load and exterior observer consume this medium. This is
//! a stationary, pressure-matched acoustic interface, not an advective mixing
//! layer, thermal equilibration or live replacement of an acoustic state.

use super::{DuctError, DuctResponse, GasState, HoleState, LossModel, Region,
    RegionalDuct, Segment, Termination, WallPin, input_impedance_load,
    input_impedance_wall, invalid};

impl RegionalDuct {
    /// Bind the homogeneous exterior independently of the internal terminal gas.
    ///
    /// This cold construction changes no geometry or stored acoustic state.
    /// Omitting it retains the historical terminal-gas exterior exactly. An
    /// explicitly supplied identical gas also takes the original uniform path.
    /// The exterior has the common static pressure but may differ in temperature
    /// or composition. Closed/ideal-open end conditions do not acquire radiation.
    ///
    /// Open/vented side holes still use the existing combined chimney/radiation
    /// law. Their gas must match the explicit exterior; a two-medium chimney is
    /// not implemented and is refused instead of using the wrong radiation.
    /// Closed holes may contain a different gas. Existing mutual-hole restrictions
    /// remain active and no second outlet load or end correction is added.
    ///
    /// # Errors
    /// Incoherent gas, unequal static pressure, or an exposed side hole whose
    /// local medium does not match this exterior.
    pub fn with_exterior_gas(mut self, gas: GasState) -> Result<Self, DuctError> {
        gas.validate_coherence().map_err(|_| invalid("exterior gas is not a coherent material state"))?;
        if gas.pressure != self.inlet_gas().pressure {
            return Err(invalid("exterior and duct need a common static pressure; mean flow is not modeled"));
        }
        for region in &self.regions {
            if region.gas == gas { continue; }
            for segment in &region.duct.segments {
                if let Segment::ToneHole { state, .. } = segment {
                    let exposed = match state {
                        HoleState::Closed => false,
                        HoleState::Open => true,
                        HoleState::Vent(sigma) => *sigma > 0.0,
                    };
                    if exposed {
                        return Err(invalid("an exposed tone hole with a different exterior gas needs an unavailable two-medium chimney law"));
                    }
                }
            }
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
        loss: LossModel, termination: Termination, wall: Option<&WallPin>)
        -> Result<DuctResponse, DuctError> {
        let exterior = self.exterior_gas();
        if *exterior == terminal.gas
            || matches!(termination, Termination::Closed | Termination::IdealOpen) {
            return input_impedance_wall(&terminal.duct, &terminal.gas,
                omega, loss, termination, wall);
        }
        let radius = terminal.duct.segments.last().expect("admitted axial outlet").outlet_radius();
        // Keep radiation with its owner, including its frequency-domain sign,
        // low-ka refusal and flanged-piston continuation. Internal propagation
        // receives that load through the same port as a measured/baked load.
        let load = fs_duct::compact_radiation_impedance(termination, exterior, radius, omega)?;
        let mut response = input_impedance_load(&terminal.duct, &terminal.gas,
            omega, loss, load, wall)?;
        response.mouth_ka = omega / exterior.sound_speed * radius;
        Ok(response)
    }
}

#[cfg(test)]
mod tests;
