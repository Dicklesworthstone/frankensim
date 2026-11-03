//! One stage's world material-air field in the original aircraft body frame.
//!
//! Hub/strip samples include carrier point motion. Blade-field differences
//! contain only material-air variation: native blade quadrature already owns
//! carrier tilt and absolute shaft rotation. Form wind differences before
//! rotation, rather than subtracting two large material-relative velocities.

use super::{Refusal, array, finite, reject, vec};
use fs_ga::{So3, Vec3};

pub(super) struct AirSample {
    pub(super) relative_velocity_mps: [f64; 3],
    pub(super) wind_world_mps: [f64; 3],
}

pub(super) struct StageIncidentAir<'a> {
    cg_body: Vec3,
    position_world: Vec3,
    velocity_world: Vec3,
    omega_body: Vec3,
    body_to_world: So3,
    world_to_body: So3,
    wind_world: &'a dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
}

impl<'a> StageIncidentAir<'a> {
    // The original mechanical owner has admitted the pose and twist.
    pub(super) fn new(
        cg_body_m: [f64; 3],
        position_world_m: [f64; 3],
        velocity_world_mps: [f64; 3],
        omega_body_rad_s: [f64; 3],
        body_to_world: So3,
        wind_world: &'a dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
    ) -> Self {
        Self {
            cg_body: vec(cg_body_m),
            position_world: vec(position_world_m),
            velocity_world: vec(velocity_world_mps),
            omega_body: vec(omega_body_rad_s),
            world_to_body: body_to_world.inverse(),
            body_to_world,
            wind_world,
        }
    }

    fn wind_at_arm(&self, arm: Vec3) -> Result<Vec3, Refusal> {
        let point_world = self.position_world
            + self.body_to_world.rotate(arm).map_err(|_| {
                reject("propulsion-state-invalid", "world sample position overflow")
            })?;
        if !finite(point_world) {
            return Err(reject(
                "propulsion-state-invalid",
                "world sample position overflow",
            ));
        }
        let wind = vec((self.wind_world)(array(point_world))?);
        if !finite(wind) {
            return Err(reject(
                "propulsion-state-invalid",
                "nonfinite material air velocity",
            ));
        }
        Ok(wind)
    }

    pub(super) fn sample(&self, point_body_m: [f64; 3]) -> Result<AirSample, Refusal> {
        let arm = vec(point_body_m) - self.cg_body;
        let wind = self.wind_at_arm(arm)?;
        // Keep the original hub/strip arithmetic, including its operation order.
        let relative = self
            .world_to_body
            .rotate(self.velocity_world - wind)
            .map_err(|_| reject("propulsion-state-invalid", "relative air velocity overflow"))?
            + self.omega_body.cross(arm);
        if !finite(relative) {
            return Err(reject(
                "propulsion-state-invalid",
                "local relative air velocity overflow",
            ));
        }
        Ok(AirSample {
            relative_velocity_mps: array(relative),
            wind_world_mps: array(wind),
        })
    }

    /// Addition to the hub-relative air velocity at a blade point, without
    /// duplicating the native blade owner's carrier or shaft point velocity.
    pub(super) fn variation(
        &self,
        point_body_m: [f64; 3],
        hub_wind_world_mps: [f64; 3],
    ) -> Result<[f64; 3], Refusal> {
        let wind = self.wind_at_arm(vec(point_body_m) - self.cg_body)?;
        let variation = self
            .world_to_body
            .rotate(vec(hub_wind_world_mps) - wind)
            .map_err(|_| reject("propulsion-state-invalid", "blade air variation overflow"))?;
        if !finite(variation) {
            return Err(reject(
                "propulsion-state-invalid",
                "blade air variation overflow",
            ));
        }
        Ok(array(variation))
    }
}
