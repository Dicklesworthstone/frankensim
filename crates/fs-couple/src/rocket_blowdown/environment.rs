//! Environment primitives shared by point and rigid rocket flight.
//!
//! Pressure and air properties come from one existing gas evaluation; wind
//! remains a velocity in the declared inertial frame. Aerodynamic speed is
//! evaluated at the actual force application point, including body rotation.
//! A supplied drag coefficient is an explicit caller model, not a hidden air
//! or rocket preset. Mach/Reynolds numbers are available to its owning law.

use fs_material::{MaterialError, fluid_transport::FluidTransport, gas::GasState};
use fs_mbd::{UnitQuaternion, Vec3};

/// Invalid frame, coordinate, constitutive state or unrepresentable load.
#[derive(Clone, Debug, PartialEq)]
pub enum FlightEnvironmentError {
    /// Invalid kinematic or force-law input.
    Invalid(&'static str),
    /// Original gas owner refused its numerical state or transport query.
    Material(MaterialError),
}
impl core::fmt::Display for FlightEnvironmentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rocket flight environment: {self:?}")
    }
}
impl std::error::Error for FlightEnvironmentError {}
impl From<MaterialError> for FlightEnvironmentError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}
type Result<T> = std::result::Result<T, FlightEnvironmentError>;
fn invalid(message: &'static str) -> FlightEnvironmentError {
    FlightEnvironmentError::Invalid(message)
}
fn norm(value: Vec3) -> f64 {
    value.x.hypot(value.y).hypot(value.z)
}

/// Complete local atmosphere and wind in a caller-identified inertial frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlightAirSample {
    frame: u64,
    gas: GasState,
    wind_world_m_s: Vec3,
}
impl FlightAirSample {
    /// Retain a coherent gas state. Source validity remains with its producer.
    pub fn new(frame: u64, gas: GasState, wind_world_m_s: Vec3) -> Result<Self> {
        gas.validate_coherence()?;
        if !wind_world_m_s.is_finite() {
            return Err(invalid("wind must be finite"));
        }
        Ok(Self { frame, gas, wind_world_m_s })
    }

    /// The same ambient pressure used for nozzle back-pressure matching [Pa].
    #[must_use]
    pub const fn pressure_pa(&self) -> f64 {
        self.gas.pressure
    }
    /// Original complete local gas snapshot.
    #[must_use]
    pub const fn gas(&self) -> &GasState {
        &self.gas
    }
    /// Ambient velocity in the declared frame [m/s].
    #[must_use]
    pub const fn wind_world_m_s(&self) -> Vec3 {
        self.wind_world_m_s
    }

    /// Air-relative flow at a point fixed in the vehicle body.
    /// The arm is measured from its center of mass. Orientation maps body
    /// coordinates to the same inertial frame. Reference length is the
    /// explicit length used by the aerodynamic coefficient model.
    pub fn at_body_point(
        &self,
        frame: u64,
        center_velocity_world_m_s: Vec3,
        orientation: UnitQuaternion,
        angular_velocity_body_rad_s: Vec3,
        arm_body_m: Vec3,
        reference_length_m: f64,
    ) -> Result<RelativeAirflow> {
        if frame != self.frame {
            return Err(invalid("air and vehicle must share the declared inertial frame"));
        }
        if !center_velocity_world_m_s.is_finite()
            || !angular_velocity_body_rad_s.is_finite()
            || !arm_body_m.is_finite()
            || !reference_length_m.is_finite()
            || reference_length_m <= 0.0
        {
            return Err(invalid("airflow needs finite kinematics and a positive reference length"));
        }
        let local_velocity =
            orientation.rotate_body_to_world(angular_velocity_body_rad_s.cross(arm_body_m));
        let relative = center_velocity_world_m_s
            .sub(self.wind_world_m_s)
            .add(local_velocity);
        let speed = norm(relative);
        let dynamic_pressure = 0.5 * self.gas.density * speed * speed;
        let mach = speed / self.gas.sound_speed;
        let reynolds = self
            .gas
            .transport_state()?
            .reynolds(speed, reference_length_m)?;
        if !relative.is_finite()
            || !speed.is_finite()
            || !dynamic_pressure.is_finite()
            || !mach.is_finite()
        {
            return Err(invalid("air-relative velocity or aerodynamic state is unrepresentable"));
        }
        Ok(RelativeAirflow {
            relative_velocity_world_m_s: relative,
            relative_velocity_body_m_s: orientation.rotate_world_to_body(relative),
            speed_m_s: speed,
            mach_number: mach,
            reynolds_number: reynolds,
            dynamic_pressure_pa: dynamic_pressure,
            orientation,
            arm_body_m,
        })
    }
}

/// Immutable coherent flow sample for the aerodynamic coefficient owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelativeAirflow {
    relative_velocity_world_m_s: Vec3,
    relative_velocity_body_m_s: Vec3,
    speed_m_s: f64,
    mach_number: f64,
    reynolds_number: f64,
    dynamic_pressure_pa: f64,
    orientation: UnitQuaternion,
    arm_body_m: Vec3,
}

/// Force and moment of a single explicitly modeled drag application point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AerodynamicLoad {
    /// Aerodynamic force in the common inertial frame [N].
    pub force_world_n: Vec3,
    /// Moment about the vehicle center of mass in body coordinates [N m].
    pub torque_body_n_m: Vec3,
    /// Force dotted with air-relative point velocity [W], nonpositive for drag.
    /// This is relative mechanical work, not an allocated vehicle heat source.
    pub relative_power_w: f64,
}
impl RelativeAirflow {
    /// Relative application-point velocity in the inertial frame [m/s].
    #[must_use]
    pub const fn relative_velocity_world_m_s(&self) -> Vec3 {
        self.relative_velocity_world_m_s
    }
    /// The same relative velocity in body coordinates [m/s].
    #[must_use]
    pub const fn relative_velocity_body_m_s(&self) -> Vec3 {
        self.relative_velocity_body_m_s
    }
    /// Relative point speed [m/s].
    #[must_use]
    pub const fn speed_m_s(&self) -> f64 {
        self.speed_m_s
    }
    /// Relative Mach number supplied to the aerodynamic model.
    #[must_use]
    pub const fn mach_number(&self) -> f64 {
        self.mach_number
    }
    /// Reynolds number for the supplied aerodynamic reference length.
    #[must_use]
    pub const fn reynolds_number(&self) -> f64 {
        self.reynolds_number
    }
    /// Dynamic pressure [Pa] from the same density and speed.
    #[must_use]
    pub const fn dynamic_pressure_pa(&self) -> f64 {
        self.dynamic_pressure_pa
    }

    /// Evaluate `-q Cd A v/|v|` and the matching application-point moment.
    /// The coefficient must come from a model valid at this Mach, Reynolds,
    /// attitude and configuration. This does not establish that validity,
    /// invent lift/stability derivatives, or predict rarefied flow.
    pub fn drag(&self, reference_area_m2: f64, drag_coefficient: f64) -> Result<AerodynamicLoad> {
        if !reference_area_m2.is_finite()
            || reference_area_m2 <= 0.0
            || !drag_coefficient.is_finite()
            || drag_coefficient < 0.0
        {
            return Err(invalid("drag needs positive finite area and nonnegative finite coefficient"));
        }
        let magnitude = self.dynamic_pressure_pa * reference_area_m2 * drag_coefficient;
        let force = if self.speed_m_s == 0.0 {
            Vec3::ZERO
        } else {
            self.relative_velocity_world_m_s
                .scale(-magnitude / self.speed_m_s)
        };
        let torque = self
            .arm_body_m
            .cross(self.orientation.rotate_world_to_body(force));
        let power = -magnitude * self.speed_m_s;
        if !force.is_finite() || !torque.is_finite() || !power.is_finite() {
            return Err(invalid("aerodynamic force, torque or work is unrepresentable"));
        }
        Ok(AerodynamicLoad {
            force_world_n: force,
            torque_body_n_m: torque,
            relative_power_w: power,
        })
    }
}

/// Explicit point-mass gravity source in a caller-declared inertial frame.
/// No terrestrial radius, reference gravity, rotation or J2 is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CentralGravity {
    frame: u64,
    center_world_m: Vec3,
    gravitational_parameter_m3_s2: f64,
}
impl CentralGravity {
    /// Bind the source position and its positive gravitational parameter G M.
    pub fn new(
        frame: u64,
        center_world_m: Vec3,
        gravitational_parameter_m3_s2: f64,
    ) -> Result<Self> {
        if !center_world_m.is_finite()
            || !gravitational_parameter_m3_s2.is_finite()
            || gravitational_parameter_m3_s2 <= 0.0
        {
            return Err(invalid("central gravity requires finite center and positive parameter"));
        }
        Ok(Self { frame, center_world_m, gravitational_parameter_m3_s2 })
    }

    /// Return acceleration [m/s²] and specific potential relative to infinity [J/kg].
    /// Repeated division avoids overflowing squared/cubed-radius intermediates.
    pub fn sample(&self, frame: u64, position_world_m: Vec3) -> Result<(Vec3, f64)> {
        if frame != self.frame {
            return Err(invalid("gravity and vehicle must share the declared inertial frame"));
        }
        let displacement = position_world_m.sub(self.center_world_m);
        let radius = norm(displacement);
        if !displacement.is_finite() || !radius.is_finite() || radius <= 0.0 {
            return Err(invalid("gravity query requires finite position away from its center"));
        }
        let potential = -self.gravitational_parameter_m3_s2 / radius;
        let acceleration = displacement.scale(1.0 / radius).scale(potential / radius);
        if !potential.is_finite() || !acceleration.is_finite() {
            return Err(invalid("central gravitational potential or acceleration is unrepresentable"));
        }
        Ok((acceleration, potential))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::gas::{GasSpec, atmosphere::StandardAtmosphere1976};

    #[test]
    fn wind_relative_drag_is_galilean_and_zero_at_co_motion() {
        let gas = StandardAtmosphere1976::new()
            .at_geometric_altitude(1000.0)
            .unwrap()
            .gas;
        let air = FlightAirSample::new(7, gas, Vec3::new(12.0, -3.0, 0.0)).unwrap();
        let velocity = Vec3::new(52.0, 0.0, 2.0);
        let sample = air
            .at_body_point(
                7, velocity, UnitQuaternion::IDENTITY, Vec3::ZERO, Vec3::ZERO, 2.0,
            )
            .unwrap();
        let load = sample.drag(0.1, 0.4).unwrap();
        assert!(load.relative_power_w < 0.0);
        let boost = Vec3::new(2000.0, -1000.0, 30.0);
        let shifted = FlightAirSample::new(7, gas, air.wind_world_m_s().add(boost))
            .unwrap()
            .at_body_point(
                7, velocity.add(boost), UnitQuaternion::IDENTITY,
                Vec3::ZERO, Vec3::ZERO, 2.0,
            )
            .unwrap();
        assert_eq!(sample, shifted);
        let rest = air
            .at_body_point(
                7, air.wind_world_m_s(), UnitQuaternion::IDENTITY,
                Vec3::ZERO, Vec3::ZERO, 2.0,
            )
            .unwrap();
        assert_eq!(rest.drag(0.1, 0.4).unwrap().force_world_n, Vec3::ZERO);
        assert_eq!(rest.mach_number(), 0.0);
        assert_eq!(rest.reynolds_number(), 0.0);
    }

    #[test]
    fn rotating_application_point_has_the_matching_dissipative_moment() {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 288.15, 101325.0).unwrap();
        let orientation =
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.7).unwrap();
        let omega = Vec3::new(0.0, 0.0, 10.0);
        let flow = FlightAirSample::new(0, gas, Vec3::ZERO)
            .unwrap()
            .at_body_point(
                0, Vec3::ZERO, orientation, omega, Vec3::new(2.0, 0.0, 0.0), 1.0,
            )
            .unwrap();
        let load = flow.drag(0.2, 0.6).unwrap();
        assert!((flow.speed_m_s() - 20.0).abs() < 1e-12);
        assert!(load.torque_body_n_m.z < 0.0);
        assert!((load.torque_body_n_m.dot(omega) - load.relative_power_w).abs() < 1e-10);
    }

    #[test]
    fn gravity_obeys_inverse_square_and_potential_gradient() {
        let gravity = CentralGravity::new(42, Vec3::new(1.0, 2.0, 3.0), 100.0).unwrap();
        let (a, potential) = gravity.sample(42, Vec3::new(3.0, 2.0, 3.0)).unwrap();
        assert_eq!(a, Vec3::new(-25.0, 0.0, 0.0));
        assert_eq!(potential, -50.0);
        let (b, _) = gravity.sample(42, Vec3::new(5.0, 2.0, 3.0)).unwrap();
        assert_eq!(b.x, a.x / 4.0);
        let h = 1e-5;
        let derivative = (gravity.sample(42, Vec3::new(3.0 + h, 2.0, 3.0)).unwrap().1
            - gravity.sample(42, Vec3::new(3.0 - h, 2.0, 3.0)).unwrap().1)
            / (2.0 * h);
        assert!((derivative + a.x).abs() < 1e-8);
        assert!(gravity.sample(0, Vec3::new(3.0, 2.0, 3.0)).is_err());
        assert!(gravity.sample(42, Vec3::new(1.0, 2.0, 3.0)).is_err());
    }
}
