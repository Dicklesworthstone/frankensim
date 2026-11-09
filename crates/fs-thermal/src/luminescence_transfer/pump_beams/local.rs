//! Externally maintained photon fields on the original finite inventory.
//!
//! These are explicit hazards per ground emitter, not transported photon
//! packets. Only actual absorbed photons supply external energy; no beam area,
//! incident/transmitted ledger or weak-depth surrogate is invented.

use super::*;

/// Local excitation hazards in inverse seconds and original photon energies
/// in joules, ordered `[ground -> singlet, ground -> triplet]`.
/// A dark/default channel may have zero photon energy. An active channel must
/// carry at least its destination state's excitation energy.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LocalPhotonPump {
    pub rates_s: [f64; 2],
    pub photon_energies_j: [f64; 2],
}

/// Actual local absorption, separate from transported and recycled photons.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LocalPumpStep {
    pub events: [f64; 2],
    pub absorbed_energy_j: [f64; 2],
    pub relaxation_heat_j: [f64; 2],
}

impl LocalPumpStep {
    pub(crate) fn from_events(
        events: [f64; 2],
        pump: LocalPhotonPump,
        state_energies: [f64; 2],
    ) -> Result<Self, SingletTripletError> {
        let mut result = Self {
            events,
            ..Self::default()
        };
        for target in 0..2 {
            result.absorbed_energy_j[target] = positive_product(
                events[target],
                pump.photon_energies_j[target],
                "local absorbed photon energy",
            )?;
            if events[target] > 0.0 {
                result.relaxation_heat_j[target] = positive_product(
                    events[target],
                    pump.photon_energies_j[target] - state_energies[target],
                    "local pump relaxation heat",
                )?;
            }
        }
        Ok(result)
    }
}

pub(super) fn validate(
    models: &[SingletTripletDecay],
    pumps: &[LocalPhotonPump],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, SingletTripletError> {
    poll(cancelled)?;
    if !pumps.is_empty() && pumps.len() != models.len() {
        return Err(SingletTripletError::Invalid("local photon pump cell count"));
    }
    let mut active = false;
    for (model, pump) in models.iter().zip(pumps) {
        poll(cancelled)?;
        for target in 0..2 {
            let rate = pump.rates_s[target];
            let energy = pump.photon_energies_j[target];
            if !nonnegative(rate)
                || !nonnegative(energy)
                || (rate > 0.0 && energy < model.excitation_energies_j()[target])
            {
                return Err(SingletTripletError::Invalid(
                    "local photon hazard or energy",
                ));
            }
            active |= rate > 0.0;
        }
    }
    Ok(active)
}

pub(super) fn at(pumps: &[LocalPhotonPump], cell: usize) -> LocalPhotonPump {
    pumps.get(cell).copied().unwrap_or_default()
}

pub(super) struct LocalTotals {
    events: [Sum; 2],
    absorbed: [Sum; 2],
    heat: [Sum; 2],
}

impl LocalTotals {
    pub(super) fn new() -> Self {
        Self {
            events: [Sum::default(); 2],
            absorbed: [Sum::default(); 2],
            heat: [Sum::default(); 2],
        }
    }

    pub(super) fn add(&mut self, step: LocalPumpStep) -> Result<(), SingletTripletError> {
        for target in 0..2 {
            self.events[target].add(step.events[target])?;
            self.absorbed[target].add(step.absorbed_energy_j[target])?;
            self.heat[target].add(step.relaxation_heat_j[target])?;
        }
        Ok(())
    }

    pub(super) fn finish(self) -> Result<LocalPumpStep, SingletTripletError> {
        Ok(LocalPumpStep {
            events: values(self.events)?,
            absorbed_energy_j: values(self.absorbed)?,
            relaxation_heat_j: values(self.heat)?,
        })
    }
}
