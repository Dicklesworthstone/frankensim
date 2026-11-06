//! Well-mixed tracer hazards from explicitly supplied carrier flows and volumes.

use fs_qty::{Volume, VolumetricFlowRate};

use super::{MAX_DECAY_NODES, NuclideOutflow, NuclideTransfer};
use crate::{FissionError, nonnegative, positive};

impl NuclideTransfer {
    /// Construct nonselective transport at `Q / V_from` for every nuclide.
    ///
    /// `flow` is an authored nonnegative volume per second; `volumes` supplies
    /// the original compartment volumes, including both referenced endpoints.
    /// The donor volume sets residence time, independently of recipient size.
    /// Zero flow retains the path with exactly zero rates for a closed valve.
    ///
    /// This describes a well-mixed tracer in fixed control volumes. It does not
    /// solve carrier continuity, change material mass, add transported enthalpy,
    /// or update the rates when geometry changes. Rebuild a rate frame explicitly
    /// when an authored flow or donor volume changes.
    pub fn from_volumetric_flow(
        from_compartment: usize,
        to_compartment: usize,
        volumes: &[Volume],
        flow: VolumetricFlowRate,
        nuclides: usize,
    ) -> Result<Self, FissionError> {
        if from_compartment == to_compartment
            || volumes
                .get(to_compartment)
                .is_none_or(|volume| !positive(volume.value()))
        {
            return Err(FissionError::InvalidInput(
                "a distinct recipient with a positive finite volume is required",
            ));
        }
        let rate = donor_rate(from_compartment, volumes, flow, nuclides)?;
        Ok(Self {
            from_compartment,
            to_compartment,
            rates_per_s: vec![rate; nuclides],
        })
    }
}

impl NuclideOutflow {
    /// Construct extraction into an untracked exterior at `Q / V_compartment`.
    ///
    /// All nuclides share the same well-mixed removal hazard. The original donor
    /// volume must be positive and finite even for zero flow. Exterior decays,
    /// carrier mass and carried sensible/latent energy remain separately owned.
    pub fn from_volumetric_flow(
        compartment: usize,
        volumes: &[Volume],
        flow: VolumetricFlowRate,
        nuclides: usize,
    ) -> Result<Self, FissionError> {
        let rate = donor_rate(compartment, volumes, flow, nuclides)?;
        Ok(Self {
            compartment,
            rates_per_s: vec![rate; nuclides],
        })
    }
}

fn donor_rate(
    compartment: usize,
    volumes: &[Volume],
    flow: VolumetricFlowRate,
    nuclides: usize,
) -> Result<f64, FissionError> {
    if nuclides == 0 || nuclides > MAX_DECAY_NODES {
        return Err(FissionError::Budget("well-mixed transport nuclide count"));
    }
    let volume = volumes
        .get(compartment)
        .filter(|volume| positive(volume.value()))
        .ok_or(FissionError::InvalidInput(
            "a positive finite donor volume is required",
        ))?
        .value();
    let flow = flow.value();
    if !nonnegative(flow) {
        return Err(FissionError::InvalidInput(
            "volumetric flow must be finite and nonnegative",
        ));
    }
    let rate = flow / volume;
    if !rate.is_finite() || (flow > 0.0 && rate == 0.0) {
        return Err(FissionError::Arithmetic);
    }
    Ok(rate)
}

#[cfg(test)]
mod tests;
