//! Prescribed material profiles on the original flight and wind lifecycle.
//!
//! Each force stage uses its CG density throughout the aircraft. The complete
//! source gas remains available at any admitted observation height; no source
//! atmosphere law or flight integration scheme is duplicated here.

use super::{SimLoop, SpatialWindSpec, decode_gas, encode_gas, invalid};
use crate::Refusal;
use crate::atmosphere::{MaterialAirProfile, Profile};
use crate::simloop::{RAIL_HEIGHT_M, ScenarioSpec};
use fs_blake3::hash_domain;

impl SimLoop {
    /// Fly through one prescribed material atmosphere with CG-sampled wind.
    /// Initialization uses the profile at the rail height, then the existing
    /// start/midpoint force stages each sample their own height. Their density
    /// drives wing, propeller, warp and lateral loads without changing the
    /// native integration, controls or wind clock.
    ///
    /// Density is frozen across the aircraft within each stage. This is not a
    /// weather, anelastic, compressible-flow or Reynolds-dependent polar solve.
    /// The initial density must satisfy the original scenario envelope.
    pub fn init_with_profile_spatial_wind(
        scenario: ScenarioSpec,
        wind: SpatialWindSpec,
        profile: MaterialAirProfile,
    ) -> Result<Self, Refusal> {
        Self::init_profile_wind(scenario, wind, profile, false)
    }

    /// The same prescribed material profile with original evolving strip/disk
    /// wind sampling. The whole aircraft shares the stage's CG density; only
    /// the incident velocity is resolved separately across the wing and rotors.
    pub fn init_with_profile_distributed_spatial_wind(
        scenario: ScenarioSpec,
        wind: SpatialWindSpec,
        profile: MaterialAirProfile,
    ) -> Result<Self, Refusal> {
        Self::init_profile_wind(scenario, wind, profile, true)
    }

    fn init_profile_wind(
        scenario: ScenarioSpec,
        wind: SpatialWindSpec,
        profile: MaterialAirProfile,
        distributed: bool,
    ) -> Result<Self, Refusal> {
        let gas = profile.gas_at_height(RAIL_HEIGHT_M)?;
        let mut sim = Self::init_material_wind(scenario, wind, gas, distributed)?;
        let mut bytes = Vec::with_capacity(161);
        bytes.extend_from_slice(sim.tick0.digest.as_bytes());
        encode_profile(&profile, &mut bytes);
        sim.tick0.digest = hash_domain("org.frankensim.wf.profile-wind-tick0.v1", &bytes).to_hex();
        bytes.clear();
        bytes.extend_from_slice(sim.run_intent_id.as_bytes());
        bytes.extend_from_slice(sim.tick0.digest.as_bytes());
        sim.run_intent_id =
            hash_domain("org.frankensim.wf.profile-wind-intent.v1", &bytes).to_hex();
        sim.digest_acc = *hash_domain(
            "org.frankensim.wf.sim-digest.v1",
            sim.tick0.digest.as_bytes(),
        )
        .as_bytes();
        let wind = sim
            .spatial_wind
            .as_mut()
            .expect("initialized spatial field");
        wind.material_gas = None;
        wind.material_profile = Some(profile);
        Ok(sim)
    }

    /// Immutable prescribed atmosphere, if explicitly bound at initialization.
    /// `scenario().rho_kg_m3` remains the initial rail density; observations and
    /// load stages use the profile at their own physical height.
    #[must_use]
    pub fn material_air_profile(&self) -> Option<&MaterialAirProfile> {
        self.spatial_wind.as_ref()?.material_profile.as_ref()
    }

    pub(in crate::simloop) fn air_density_at_height(&self, height_m: f64) -> Result<f64, Refusal> {
        match self.material_air_profile() {
            Some(profile) => profile.gas_at_height(height_m).map(|gas| gas.density),
            None => Ok(self.spec.rho_kg_m3),
        }
    }
}

// V6 stores the model definition, not a gas snapshot frozen at the checkpoint
// altitude. The profile identity binds the source law and geometric datum.
pub(in crate::simloop) fn encode_profile(profile: &MaterialAirProfile, bytes: &mut Vec<u8>) {
    match &profile.profile {
        Profile::Uniform(gas) => {
            bytes.push(0);
            encode_gas(gas, bytes);
        }
        Profile::Standard1976 {
            ground_altitude_m, ..
        } => {
            bytes.push(1);
            bytes.extend_from_slice(&ground_altitude_m.to_bits().to_le_bytes());
        }
    }
}

pub(in crate::simloop) fn decode_profile(
    bytes: &[u8],
) -> Result<(MaterialAirProfile, usize), Refusal> {
    match bytes.first() {
        Some(0) => {
            let gas = decode_gas(
                bytes
                    .get(1..97)
                    .ok_or_else(|| invalid("truncated uniform atmosphere profile"))?,
            )?;
            Ok((MaterialAirProfile::uniform(gas)?, 97))
        }
        Some(1) => {
            let ground = f64::from_le_bytes(
                bytes
                    .get(1..9)
                    .ok_or_else(|| invalid("truncated hydrostatic atmosphere profile"))?
                    .try_into()
                    .expect("eight-byte ground altitude"),
            );
            Ok((MaterialAirProfile::standard1976(ground)?, 9))
        }
        _ => Err(invalid("unknown material atmosphere profile")),
    }
}

#[cfg(test)]
mod tests;
