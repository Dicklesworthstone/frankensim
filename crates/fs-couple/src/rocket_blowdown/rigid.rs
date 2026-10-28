//! One source-bearing reservoir drives the original rigid rocket owner.
//!
//! The gas owner is solved once per interval. Its accepted expelled mass and
//! scalar axial impulse determine the body-directed effective exhaust velocity;
//! the mechanics owner resolves attitude and inertia on its own bounded substeps.
//! Rotating thrust has a world impulse whose norm can be smaller than the axial
//! impulse. Only their scalar axial quantities may be compared for conservation.
//! Pressure thrust enters both the net thrust and its lever-arm torque once.
//!
//! Gas thermodynamics is comoving and spatially mixed. The native rigid model
//! retains its fixed-COM/common-principal-axis inertia law and explicitly supplied
//! diagonal rotational exhaust flux. This adapter does not infer internal-flow
//! momentum, general outlet tensors, slosh, COM motion or world exhaust energy.

use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::blowdown::{
    GasBlowdownError, GasBlowdownStep, RocketBlowdownSolve,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::{
    RocketNozzle, RocketNozzleFlow,
};
use fs_mbd::rocket::RocketState;
use fs_mbd::rocket::rigid::{
    RigidRocketError, RigidRocketForcing, RigidRocketInertia, RigidRocketSolve, RigidRocketState,
    RigidRocketStep,
};
use fs_mbd::{UnitQuaternion, Vec3};

/// Forcing in one caller-declared inertial frame and the vehicle body frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasRigidRocketEnvironment {
    /// Identity of the common world frame; no coordinate transform is inferred.
    pub inertial_frame: u64,
    /// Explicit nonnegative back pressure [Pa], zero for vacuum.
    pub ambient_pressure_pa: f64,
    /// Unit body-frame vector pointing along total nozzle thrust.
    pub thrust_direction_body: Vec3,
    /// Effective thrust application point relative to the shared COM [m].
    pub thrust_application_arm_body_m: Vec3,
    /// Declared diagonal outflow co-rotation moment per mass [m²].
    /// Independent of the thrust arm; arbitrary off-axis outlet tensors are unsupported.
    pub co_rotating_exhaust_specific_inertia_m2: Vec3,
    /// Additional exported angular momentum per mass [m²/s], excluding axial jet torque.
    pub exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3,
    /// Gravitational acceleration in the matching world frame [m/s²].
    pub gravity_world_m_s2: Vec3,
    /// Other world-frame force [N], excluding gravity and net nozzle thrust.
    pub external_force_world_n: Vec3,
    /// Other body-frame torque about the COM [N m].
    pub external_torque_body_n_m: Vec3,
}

/// Separate controls for original gas and rigid-body owners.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasRigidRocketSolve {
    /// Source, reservoir and finite-time nozzle integration controls.
    pub gas: RocketBlowdownSolve,
    /// Attitude/direction resolution and native substep budget.
    pub mechanics: RigidRocketSolve,
}

/// Loads during an explicitly closed-valve interval; no nozzle solve is implied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasRigidRocketCoastEnvironment {
    /// Identity of the common inertial frame.
    pub inertial_frame: u64,
    /// World gravitational acceleration [m/s²].
    pub gravity_world_m_s2: Vec3,
    /// Other world force [N], excluding gravity and any engine thrust.
    pub external_force_world_n: Vec3,
    /// Body-frame torque about the COM [N m].
    pub external_torque_body_n_m: Vec3,
}

/// An immutable gas reservoir and vehicle whose propellant comes only from it.
#[derive(Clone, Debug, PartialEq)]
pub struct GasRigidRocketState {
    reservoir: ClosedGasParcel,
    vehicle: RigidRocketState,
    inertial_frame: u64,
}

/// Refusal without publishing either new physical owner.
#[derive(Clone, Debug, PartialEq)]
pub enum GasRigidRocketError {
    /// A supplied frame or direction is invalid.
    Invalid(&'static str),
    /// The source/reservoir/nozzle owner refused the physical interval.
    Gas(GasBlowdownError),
    /// The original native rigid rocket owner refused.
    Mechanics(RigidRocketError),
    /// Actual consumed or remaining inventories differ beyond the gas mass tolerance.
    InventoryBalance {
        /// Vehicle remaining propellant minus remaining tank gas [kg].
        remaining_residual_kg: f64,
        /// Accepted native consumed mass minus actual expelled gas [kg].
        consumed_residual_kg: f64,
    },
    /// Actual scalar mechanical axial impulse differs from the accepted gas impulse.
    AxialImpulseBalance {
        /// Mechanical minus gas axial impulse [N s].
        residual_n_s: f64,
    },
    /// Cancellation at the coupling boundary.
    Cancelled,
}
impl core::fmt::Display for GasRigidRocketError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-gas rigid rocket coupling: {self:?}")
    }
}
impl std::error::Error for GasRigidRocketError {}
impl From<GasBlowdownError> for GasRigidRocketError {
    fn from(value: GasBlowdownError) -> Self {
        Self::Gas(value)
    }
}
impl From<RigidRocketError> for GasRigidRocketError {
    fn from(value: RigidRocketError) -> Self {
        Self::Mechanics(value)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasRigidRocketError> {
    if cancelled() {
        Err(GasRigidRocketError::Cancelled)
    } else {
        Ok(())
    }
}

/// Original gas/nozzle and native mechanics reports for the same interval.
#[derive(Clone, Debug, PartialEq)]
pub struct GasRigidRocketStep {
    /// Shared physical duration [s].
    pub duration_s: f64,
    /// Actual finite reservoir mass, species and stagnation-enthalpy accounting.
    pub gas: GasBlowdownStep,
    /// Original source-resolved midpoint nozzle flow and regime.
    pub midpoint_nozzle: RocketNozzleFlow,
    /// Original native rigid trajectory, directed impulse and spin accounting.
    pub mechanics: RigidRocketStep,
    /// Remaining vehicle propellant minus remaining tank mass [kg].
    pub remaining_mass_residual_kg: f64,
    /// Accepted native consumed mass minus actual expelled gas [kg].
    pub consumed_mass_residual_kg: f64,
    /// Mechanical scalar axial impulse minus accepted nozzle impulse [N s].
    pub axial_impulse_residual_n_s: f64,
}

/// Complete candidate gas and rigid endpoints, published atomically.
#[derive(Clone, Debug, PartialEq)]
pub struct GasRigidRocketCandidate {
    endpoint: GasRigidRocketState,
    step: GasRigidRocketStep,
}

/// Complete closed-valve candidate; reservoir mass and thermodynamic state are unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct GasRigidRocketCoastCandidate {
    endpoint: GasRigidRocketState,
    mechanics: RigidRocketStep,
}
impl GasRigidRocketCoastCandidate {
    /// Complete next state with the identical source-bearing reservoir.
    #[must_use]
    pub const fn endpoint(&self) -> &GasRigidRocketState {
        &self.endpoint
    }
    /// Original native rigid coast report.
    #[must_use]
    pub const fn mechanics(&self) -> &RigidRocketStep {
        &self.mechanics
    }
    /// Accept coast and its native report without a gas query or another step.
    #[must_use]
    pub fn into_parts(self) -> (GasRigidRocketState, RigidRocketStep) {
        (self.endpoint, self.mechanics)
    }
}
impl GasRigidRocketCandidate {
    /// Joint candidate endpoint.
    #[must_use]
    pub const fn endpoint(&self) -> &GasRigidRocketState {
        &self.endpoint
    }
    /// Original owner reports and actual coupling residuals.
    #[must_use]
    pub const fn step(&self) -> &GasRigidRocketStep {
        &self.step
    }
    /// Accept both physical owners without another solve or inventory repair.
    #[must_use]
    pub fn into_parts(self) -> (GasRigidRocketState, GasRigidRocketStep) {
        (self.endpoint, self.step)
    }
}

impl GasRigidRocketState {
    /// Assemble native rigid motion with all propellant in this source reservoir.
    /// Dry mass includes the tank/container. Propellant inventory is derived.
    pub fn new(
        reservoir: ClosedGasParcel,
        dry_mass_kg: f64,
        position_m: Vec3,
        velocity_m_s: Vec3,
        orientation: UnitQuaternion,
        angular_momentum_body_kg_m2_s: Vec3,
        inertia: RigidRocketInertia,
        inertial_frame: u64,
    ) -> Result<Self, GasRigidRocketError> {
        let translation = RocketState {
            position_m,
            velocity_m_s,
            dry_mass_kg,
            propellant_mass_kg: reservoir.mass_kg(),
        };
        let vehicle = RigidRocketState::new(
            translation,
            orientation,
            angular_momentum_body_kg_m2_s,
            inertia,
        )?;
        Ok(Self {
            reservoir,
            vehicle,
            inertial_frame,
        })
    }
    /// Original evolving source-bearing gas state.
    #[must_use]
    pub const fn reservoir(&self) -> &ClosedGasParcel {
        &self.reservoir
    }
    /// Original native rigid rocket state.
    #[must_use]
    pub const fn vehicle(&self) -> &RigidRocketState {
        &self.vehicle
    }
    /// Identity of the common world frame, not proof that it is inertial.
    #[must_use]
    pub const fn inertial_frame(&self) -> u64 {
        self.inertial_frame
    }

    /// Advance an explicitly closed valve under the supplied external loads.
    ///
    /// Engine mass flow, thrust and rotational exhaust flux are exactly zero.
    /// The insulated fixed-volume reservoir stays unchanged. This is an explicit
    /// control mode, never a fallback for a failed source/nozzle solve. Existing
    /// admitted inventory roundoff is retained, not hidden by resetting an owner.
    pub fn coast(
        &self,
        environment: GasRigidRocketCoastEnvironment,
        duration_s: f64,
        controls: RigidRocketSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasRigidRocketCoastCandidate, GasRigidRocketError> {
        poll(&mut cancelled)?;
        if environment.inertial_frame != self.inertial_frame {
            return Err(GasRigidRocketError::Invalid(
                "rigid coast state and loads must share the declared inertial frame",
            ));
        }
        let mechanics = self.vehicle.advance(
            RigidRocketForcing {
                mass_flow_kg_s: 0.0,
                effective_exhaust_velocity_body_m_s: Vec3::ZERO,
                thrust_application_arm_body_m: Vec3::ZERO,
                co_rotating_exhaust_specific_inertia_m2: Vec3::ZERO,
                exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3::ZERO,
                gravity_world_m_s2: environment.gravity_world_m_s2,
                external_force_world_n: environment.external_force_world_n,
                external_torque_body_n_m: environment.external_torque_body_n_m,
            },
            duration_s,
            controls,
            &mut cancelled,
        )?;
        if mechanics.burned_propellant_mass_kg != 0.0
            || mechanics.axial_thrust_impulse_n_s != 0.0
            || mechanics.state_after.translation().propellant_mass_kg
                != self.vehicle.translation().propellant_mass_kg
            || mechanics.state_after.translation().dry_mass_kg
                != self.vehicle.translation().dry_mass_kg
        {
            return Err(GasRigidRocketError::Invalid(
                "closed-valve native coast changed an engine inventory or impulse",
            ));
        }
        poll(&mut cancelled)?;
        let endpoint = Self {
            reservoir: self.reservoir.clone(),
            vehicle: mechanics.state_after,
            inertial_frame: self.inertial_frame,
        };
        Ok(GasRigidRocketCoastCandidate {
            endpoint,
            mechanics,
        })
    }

    /// Advance the gas once and consume its accepted extensive quantities.
    ///
    /// Mechanics substeps refine rotating body forcing, not the gas quadrature.
    /// Refine the shared interval for changing nozzle conditions/environment.
    /// Native inventory rounding remains visible and is admitted only within
    /// the source reservoir mass tolerance; no owner is silently resynchronized.
    pub fn advance(
        &self,
        environment: GasRigidRocketEnvironment,
        nozzle: RocketNozzle,
        duration_s: f64,
        controls: GasRigidRocketSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasRigidRocketCandidate, GasRigidRocketError> {
        poll(&mut cancelled)?;
        if environment.inertial_frame != self.inertial_frame {
            return Err(GasRigidRocketError::Invalid(
                "rigid gas state and forcing must share the declared inertial frame",
            ));
        }
        let direction = environment.thrust_direction_body;
        let norm_squared = direction.norm_squared();
        if !direction.is_finite()
            || !norm_squared.is_finite()
            || (norm_squared - 1.0).abs() > 64.0 * f64::EPSILON
        {
            return Err(GasRigidRocketError::Invalid(
                "body thrust direction must be a finite unit vector",
            ));
        }
        let mut forcing = RigidRocketForcing {
            mass_flow_kg_s: 0.0,
            effective_exhaust_velocity_body_m_s: Vec3::ZERO,
            thrust_application_arm_body_m: environment.thrust_application_arm_body_m,
            co_rotating_exhaust_specific_inertia_m2: environment
                .co_rotating_exhaust_specific_inertia_m2,
            exhaust_swirl_specific_angular_momentum_body_m2_s: environment
                .exhaust_swirl_specific_angular_momentum_body_m2_s,
            gravity_world_m_s2: environment.gravity_world_m_s2,
            external_force_world_n: environment.external_force_world_n,
            external_torque_body_n_m: environment.external_torque_body_n_m,
        };
        // Validate native controls and finite forcing before the expensive gas solve.
        self.vehicle
            .advance(forcing, 0.0, controls.mechanics, &mut cancelled)?;
        let gas = self.reservoir.rocket_blowdown(
            environment.ambient_pressure_pa,
            nozzle,
            duration_s,
            controls.gas,
            &mut cancelled,
        )?;
        let expelled = gas.reservoir().step().expelled_mass_kg;
        let axial_impulse = gas.axial_impulse_n_s;
        let rate = expelled / duration_s;
        let exhaust = axial_impulse / expelled;
        if !rate.is_finite() || rate <= 0.0 || !exhaust.is_finite() || exhaust <= 0.0 {
            return Err(GasRigidRocketError::Invalid(
                "accepted expelled mass and axial impulse must define finite positive forcing",
            ));
        }
        forcing.mass_flow_kg_s = rate;
        forcing.effective_exhaust_velocity_body_m_s = direction.scale(exhaust);
        let mechanics =
            self.vehicle
                .advance(forcing, duration_s, controls.mechanics, &mut cancelled)?;
        let remaining_residual = mechanics.state_after.translation().propellant_mass_kg
            - gas.reservoir().endpoint().mass_kg();
        let consumed_residual = mechanics.burned_propellant_mass_kg - expelled;
        let tolerance = controls.gas.reservoir.mass_tolerance_kg;
        if !remaining_residual.is_finite()
            || !consumed_residual.is_finite()
            || remaining_residual.abs() > tolerance
            || consumed_residual.abs() > tolerance
            || mechanics.propellant_balance_residual_kg.abs() > tolerance
            || mechanics.coast_duration_s != 0.0
            || mechanics.state_after.translation().dry_mass_kg
                != self.vehicle.translation().dry_mass_kg
        {
            return Err(GasRigidRocketError::InventoryBalance {
                remaining_residual_kg: remaining_residual,
                consumed_residual_kg: consumed_residual,
            });
        }
        // World impulse is a vector integral over attitude. Its norm is NOT
        // the conserved axial scalar and must not replace this comparison.
        let impulse_residual = mechanics.axial_thrust_impulse_n_s - axial_impulse;
        if !impulse_residual.is_finite()
            || impulse_residual.abs() > 128.0 * f64::EPSILON * axial_impulse
        {
            return Err(GasRigidRocketError::AxialImpulseBalance {
                residual_n_s: impulse_residual,
            });
        }
        poll(&mut cancelled)?;
        let (reservoir, gas, midpoint_nozzle) = gas.into_parts();
        let endpoint = Self {
            reservoir,
            vehicle: mechanics.state_after,
            inertial_frame: self.inertial_frame,
        };
        Ok(GasRigidRocketCandidate {
            endpoint,
            step: GasRigidRocketStep {
                duration_s,
                gas,
                midpoint_nozzle,
                mechanics,
                remaining_mass_residual_kg: remaining_residual,
                consumed_mass_residual_kg: consumed_residual,
                axial_impulse_residual_n_s: impulse_residual,
            },
        })
    }
}

#[cfg(test)]
#[path = "rigid/tests.rs"]
mod tests;
