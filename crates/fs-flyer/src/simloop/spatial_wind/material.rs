//! Shared material gas on the existing prescribed spatial-wind lifecycle.
//!
//! The complete numerical gas state is retained without re-evaluating its law.
//! Its density drives the original trim, wing, propeller and motion owners.
//! The current polar model has no Reynolds/temperature/compressibility axis;
//! carrying those properties does not invent their effect on the polar.

use super::{FlowSample, SimLoop, SpatialWindSpec};
use crate::Refusal;
use crate::simloop::ScenarioSpec;
use fs_blake3::hash_domain;
use fs_material::gas::GasState;
use fs_math::det;

/// The accepted wind and its complete, homogeneous material gas at one tick.
/// Components and gradient coordinates are flight-line forward/right/up.
/// Numerical gas identity is retained; source receipts stay with the caller's
/// material model, and this view grants no composition or calibration authority.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialAirSample {
    /// The exact original wind sample, without changing its random stream.
    pub wind: FlowSample,
    /// Original admitted thermodynamics and transport, with no air substitution.
    pub gas: GasState,
}

impl MaterialAirSample {
    fn relative_speed_squared(&self, observer_velocity_mps: [f64; 3]) -> Result<f64, Refusal> {
        if !observer_velocity_mps.iter().all(|v| v.is_finite()) {
            return Err(invalid(
                "air-relative diagnostics need finite observer velocity",
            ));
        }
        let speed2: f64 = self
            .wind
            .u
            .iter()
            .zip(observer_velocity_mps)
            .map(|(wind, observer)| (wind - observer) * (wind - observer))
            .sum();
        if !speed2.is_finite() {
            return Err(invalid("air-relative speed is unrepresentable"));
        }
        Ok(speed2)
    }

    /// Dynamic pressure relative to an explicit observer velocity [Pa].
    /// A moving aircraft uses its ground velocity, not the ambient wind speed.
    pub fn dynamic_pressure_pa(&self, observer_velocity_mps: [f64; 3]) -> Result<f64, Refusal> {
        let q = 0.5 * self.gas.density * self.relative_speed_squared(observer_velocity_mps)?;
        if !q.is_finite() {
            return Err(invalid("dynamic pressure is unrepresentable"));
        }
        Ok(q)
    }

    /// Reynolds number from the same local relative velocity, density and mu.
    /// This is a diagnostic; the existing reference polar is Reynolds-independent.
    pub fn reynolds(&self, chord_m: f64, observer_velocity_mps: [f64; 3]) -> Result<f64, Refusal> {
        if !chord_m.is_finite() || chord_m <= 0.0 {
            return Err(invalid("Reynolds number needs a positive finite chord"));
        }
        let re = self.gas.density
            * det::sqrt(self.relative_speed_squared(observer_velocity_mps)?)
            * chord_m
            / self.gas.dynamic_viscosity;
        if !re.is_finite() {
            return Err(invalid("Reynolds number is unrepresentable"));
        }
        Ok(re)
    }
}

fn invalid(message: &str) -> Refusal {
    Refusal {
        code: "material-atmosphere-invalid",
        message: message.into(),
        ranked_repairs: vec!["evaluate one coherent GasState through its material owner".into()],
    }
}

impl SimLoop {
    /// Initialize CG-sampled wind with one complete material gas. Its density
    /// replaces the scalar scenario density BEFORE native trim and initialization.
    /// The accepted scenario is available from `scenario()` for exact restoration.
    /// Gas properties are homogeneous and prescribed for the whole run; heat,
    /// composition evolution and Reynolds-dependent polar corrections are absent.
    /// The existing density envelope (0.5, 2.0) kg/m3 still applies.
    pub fn init_with_material_spatial_wind(
        scenario: ScenarioSpec,
        wind: SpatialWindSpec,
        gas: GasState,
    ) -> Result<Self, Refusal> {
        Self::init_material_wind(scenario, wind, gas, false)
    }

    /// The same material state with the original distributed wing/canard/disk
    /// wind sampling. Trim, all force paths and rotor torque share its density.
    /// All approximation boundaries of the distributed spatial tier remain.
    pub fn init_with_material_distributed_spatial_wind(
        scenario: ScenarioSpec,
        wind: SpatialWindSpec,
        gas: GasState,
    ) -> Result<Self, Refusal> {
        Self::init_material_wind(scenario, wind, gas, true)
    }

    fn init_material_wind(
        mut scenario: ScenarioSpec,
        wind: SpatialWindSpec,
        gas: GasState,
        distributed: bool,
    ) -> Result<Self, Refusal> {
        gas.validate_coherence()
            .map_err(|e| invalid(&e.to_string()))?;
        scenario.rho_kg_m3 = gas.density;
        let mut sim = if distributed {
            Self::init_with_distributed_spatial_wind(scenario, wind)?
        } else {
            Self::init_with_spatial_wind(scenario, wind)?
        };
        let mut bytes = Vec::with_capacity(160);
        bytes.extend_from_slice(sim.tick0.digest.as_bytes());
        encode_gas(&gas, &mut bytes);
        sim.tick0.digest = hash_domain("org.frankensim.wf.material-wind-tick0.v1", &bytes).to_hex();
        bytes.clear();
        bytes.extend_from_slice(sim.run_intent_id.as_bytes());
        bytes.extend_from_slice(sim.tick0.digest.as_bytes());
        sim.run_intent_id =
            hash_domain("org.frankensim.wf.material-wind-intent.v1", &bytes).to_hex();
        sim.digest_acc = *hash_domain(
            "org.frankensim.wf.sim-digest.v1",
            sim.tick0.digest.as_bytes(),
        )
        .as_bytes();
        sim.spatial_wind
            .as_mut()
            .expect("initialized spatial field")
            .material_gas = Some(gas);
        Ok(sim)
    }

    /// Original complete gas, if this lifecycle explicitly bound a material.
    /// Historical scalar-density scenarios have no inferred gas composition.
    #[must_use]
    pub fn material_gas(&self) -> Option<&GasState> {
        self.spatial_wind.as_ref()?.material_gas.as_ref()
    }

    /// Read wind, thermodynamics and transport from the same accepted state.
    /// This observation neither advances physics nor consumes randomness.
    pub fn sample_material_air(&self, position_m: [f64; 3]) -> Result<MaterialAirSample, Refusal> {
        let gas = *self.material_gas().ok_or_else(|| {
            invalid("this run has no material-bound gas; initialize with material spatial wind")
        })?;
        Ok(MaterialAirSample {
            wind: self.sample_spatial_wind(position_m)?,
            gas,
        })
    }
}

// A numerical snapshot, not serialization of source authority. Keep all twelve
// fields in the shared GasState order so no observer loses transport/caloric data.
pub(crate) fn encode_gas(gas: &GasState, bytes: &mut Vec<u8>) {
    for value in [
        gas.temperature,
        gas.pressure,
        gas.density,
        gas.sound_speed,
        gas.dynamic_viscosity,
        gas.thermal_conductivity,
        gas.gamma,
        gas.specific_gas_constant,
        gas.specific_heat_cp,
        gas.prandtl,
        gas.characteristic_impedance,
        gas.water_mole_fraction,
    ] {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
}

pub(crate) fn decode_gas(bytes: &[u8]) -> Result<GasState, Refusal> {
    if bytes.len() != 96 {
        return Err(invalid("truncated material gas snapshot"));
    }
    let mut fields = [0.0; 12];
    for (field, chunk) in fields.iter_mut().zip(bytes.chunks_exact(8)) {
        *field = f64::from_le_bytes(chunk.try_into().expect("eight-byte gas field"));
    }
    let [
        temperature,
        pressure,
        density,
        sound_speed,
        dynamic_viscosity,
        thermal_conductivity,
        gamma,
        specific_gas_constant,
        specific_heat_cp,
        prandtl,
        characteristic_impedance,
        water_mole_fraction,
    ] = fields;
    let gas = GasState {
        temperature,
        pressure,
        density,
        sound_speed,
        dynamic_viscosity,
        thermal_conductivity,
        gamma,
        specific_gas_constant,
        specific_heat_cp,
        prandtl,
        characteristic_impedance,
        water_mole_fraction,
    };
    gas.validate_coherence()
        .map_err(|e| invalid(&e.to_string()))?;
    Ok(gas)
}

#[cfg(test)]
mod tests;
