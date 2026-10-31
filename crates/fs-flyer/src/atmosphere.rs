//! Material gas and hydrostatic pressure on the original atmospheric wind.
//!
//! Composition happens at L4: the L2 wind owner does not depend on L3 materials.
//! The standard profile remains the existing material owner's USSA-1976 model;
//! no pressure, density, viscosity or sound-speed formula is duplicated here.
//! The optional force path samples this same state before the native aircraft
//! and propeller solve. It does not alter the historical homogeneous SimLoop.

use crate::Refusal;
use crate::aircraft::{ForceBuildup, OpenLoopDesign};
use crate::simloop::MaterialAirSample;
use fs_atmo::evolving::EvolvingAtmosphere;
use fs_atmo::{Atmosphere, FlowSample};
use fs_material::gas::GasState;
use fs_material::gas::atmosphere::StandardAtmosphere1976;
use fs_math::det;

#[derive(Clone, Copy, Debug, PartialEq)]
enum Profile {
    Uniform(GasState),
    Standard1976 {
        model: StandardAtmosphere1976,
        ground_altitude_m: f64,
    },
}

/// One immutable gas model, independent of the wind's legacy constant air data.
/// Standard heights are geometric sea-level altitudes; wind heights remain
/// distance above the aerodynamic ground plane. Their offset is explicit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialAirProfile {
    profile: Profile,
}

fn invalid(message: impl Into<String>) -> Refusal {
    Refusal {
        code: "material-atmosphere-invalid",
        message: message.into(),
        ranked_repairs: vec![
            "use a coherent material gas and a height within its declared profile".into(),
        ],
    }
}

fn wind_refusal(error: fs_atmo::Refusal) -> Refusal {
    Refusal {
        code: error.code,
        message: error.message,
        ranked_repairs: error.ranked_repairs,
    }
}

impl MaterialAirProfile {
    /// Retain any coherent, homogeneous material gas, including a non-air
    /// mixture. This does not infer composition authority from numerical data;
    /// source receipts remain with the caller's original material owner.
    pub fn uniform(gas: GasState) -> Result<Self, Refusal> {
        gas.validate_coherence()
            .map_err(|error| invalid(error.to_string()))?;
        Ok(Self {
            profile: Profile::Uniform(gas),
        })
    }

    /// Use the original dry-air USSA-1976 hydrostatic profile. The material
    /// source admits geometric sea-level altitude in [0, 80000] m; the given
    /// ground elevation and every later query must remain in that domain.
    /// This is a reference atmosphere, not observed local weather or humidity.
    pub fn standard1976(ground_altitude_m: f64) -> Result<Self, Refusal> {
        let model = StandardAtmosphere1976::new();
        model
            .at_geometric_altitude(ground_altitude_m)
            .map_err(|error| invalid(error.to_string()))?;
        Ok(Self {
            profile: Profile::Standard1976 {
                model,
                ground_altitude_m,
            },
        })
    }

    /// Complete coherent thermodynamics and transport at height above ground.
    /// A query neither changes wind nor advances any material state.
    pub fn gas_at_height(&self, height_m: f64) -> Result<GasState, Refusal> {
        if !height_m.is_finite() || height_m < 0.0 {
            return Err(invalid(
                "material air needs finite nonnegative height above ground",
            ));
        }
        match &self.profile {
            Profile::Uniform(gas) => Ok(*gas),
            Profile::Standard1976 {
                model,
                ground_altitude_m,
            } => model
                .at_geometric_altitude(*ground_altitude_m + height_m)
                .map(|state| state.gas)
                .map_err(|error| invalid(error.to_string())),
        }
    }

    /// Sample frozen atmospheric wind and this gas at the same physical point.
    /// Positions and output vectors are flight-line forward/right/up; forward
    /// is opposite the atmospheric owner's along-wind coordinate. Both vector
    /// and gradient bases are transformed, preserving analytic derivatives.
    pub fn sample(
        &self,
        wind: &Atmosphere,
        [x, y, h]: [f64; 3],
        tick: u64,
    ) -> Result<MaterialAirSample, Refusal> {
        let gas = self.gas_at_height(h)?;
        let flow = wind.sample_flow(-x, y, h, tick).map_err(wind_refusal)?;
        Ok(MaterialAirSample {
            wind: flight_frame(flow),
            gas,
        })
    }

    /// The same material sampling at an evolving atmosphere's current tick.
    /// No OU advance, random draw or checkpoint mutation occurs here.
    pub fn sample_evolving(
        &self,
        wind: &EvolvingAtmosphere,
        [x, y, h]: [f64; 3],
    ) -> Result<MaterialAirSample, Refusal> {
        let gas = self.gas_at_height(h)?;
        let flow = wind.sample_flow(-x, y, h).map_err(wind_refusal)?;
        Ok(MaterialAirSample {
            wind: flight_frame(flow),
            gas,
        })
    }
}

fn flight_frame(mut flow: FlowSample) -> FlowSample {
    let signs = [-1.0, 1.0, 1.0];
    for i in 0..3 {
        flow.u[i] *= signs[i];
        for j in 0..3 {
            flow.grad[i][j] *= signs[i] * signs[j];
        }
    }
    flow
}

/// Explicit state for the existing reduced longitudinal aircraft load model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtmosphericFlightState {
    /// Flight-line forward/right/up position [m], with height above ground.
    pub position_m: [f64; 3],
    /// Inertial forward/up ground velocity [m/s]. Wind is subtracted once.
    pub ground_velocity_mps: [f64; 2],
    /// Nose-up pitch relative to the horizontal [rad].
    pub pitch_rad: f64,
    /// Canard deflection [rad].
    pub canard_rad: f64,
    /// Propeller angular speed [rad/s].
    pub propeller_rad_s: f64,
    /// Body pitch angular speed [rad/s].
    pub pitch_rate_rad_s: f64,
}

/// Actual original gas/wind sample beside the loads that consumed its density.
#[derive(Clone, Debug)]
pub struct AtmosphericLoads {
    /// Original material state and wind sampled at the requested CG position.
    pub air: MaterialAirSample,
    /// Longitudinal air-relative speed actually supplied to the native loads.
    pub airspeed_mps: f64,
    /// Body air incidence actually supplied to the native loads [rad].
    pub angle_of_attack_rad: f64,
    /// Original aircraft/propeller result, including its level-flight
    /// `theta = alpha` gravity convention. An inertial integrator removes that
    /// contribution and applies its declared world gravity, as SimLoop does.
    pub loads: ForceBuildup,
}

impl OpenLoopDesign {
    /// Evaluate original wing, drag, propeller and engine loads using the same
    /// profile pressure/density and locally sampled relative wind. Sampling is
    /// at the CG and density is frozen across the aircraft for this load solve.
    /// Crosswind is retained for observation but the original longitudinal model
    /// has no sideslip response. Reynolds/Mach diagnostics do not change the
    /// original polar; no compressible-flow or high-altitude validation is implied.
    /// The returned load frame and level-flight gravity convention are those
    /// of `force_buildup`; this method does not itself integrate flight motion.
    pub fn force_buildup_in_atmosphere(
        &self,
        profile: &MaterialAirProfile,
        wind: &Atmosphere,
        state: AtmosphericFlightState,
        tick: u64,
    ) -> Result<AtmosphericLoads, Refusal> {
        self.material_air_loads(profile.sample(wind, state.position_m, tick)?, state)
    }

    /// Evaluate the same load path at the evolving atmosphere's accepted clock.
    /// The atmosphere is only borrowed; this call cannot advance its turbulence.
    pub fn force_buildup_in_evolving_atmosphere(
        &self,
        profile: &MaterialAirProfile,
        wind: &EvolvingAtmosphere,
        state: AtmosphericFlightState,
    ) -> Result<AtmosphericLoads, Refusal> {
        self.material_air_loads(profile.sample_evolving(wind, state.position_m)?, state)
    }

    fn material_air_loads(
        &self,
        air: MaterialAirSample,
        state: AtmosphericFlightState,
    ) -> Result<AtmosphericLoads, Refusal> {
        if !state.ground_velocity_mps.iter().all(|v| v.is_finite()) || !state.pitch_rad.is_finite()
        {
            return Err(invalid(
                "atmospheric loads need finite ground velocity and pitch",
            ));
        }
        let [u, w] = crate::simloop::motion::body_air_velocity(
            state.ground_velocity_mps,
            state.pitch_rad,
            [air.wind.u[0], air.wind.u[2]],
        );
        let speed = det::sqrt(u * u + w * w);
        let alpha = det::atan2(w, u);
        let loads = self.force_buildup(
            speed,
            alpha,
            state.canard_rad,
            state.propeller_rad_s,
            state.pitch_rate_rad_s,
            air.gas.density,
        )?;
        Ok(AtmosphericLoads {
            air,
            airspeed_mps: speed,
            angle_of_attack_rad: alpha,
            loads,
        })
    }
}

#[cfg(test)]
mod tests;
