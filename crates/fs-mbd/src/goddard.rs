//! Source-bounded kinematics for Goddard's 1914 apparatus plus the retained
//! adjacent liquid-rocket teaching kernel used by the later catalogue record.

use crate::{
    DynamicsError, Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState,
    UnitQuaternion, Vec3, Wrench,
};

/// Claim 2's stated lower bound for tapered-tube length divided by its
/// longest diameter.
pub const GODDARD_CLAIM_2_MIN_TUBE_LENGTH_RATIO: f64 = 3.0;

/// Lowest tube ratio admitted by the break-the-claim teaching surface.
pub const GODDARD_MIN_TUBE_LENGTH_RATIO: f64 = 1.0;

/// Highest tube ratio admitted by the bounded browser surface.
pub const GODDARD_MAX_TUBE_LENGTH_RATIO: f64 = 12.0;

/// Maximum elapsed time admitted by one browser pose query.
pub const GODDARD_MAX_ELAPSED_SECONDS: f64 = 600.0;

/// Maximum primary-rocket display spin admitted by the browser surface.
pub const GODDARD_MAX_PRIMARY_SPIN_RPM: f64 = 1_200.0;

/// Maximum gyroscope display spin admitted by the browser surface.
pub const GODDARD_MAX_GYRO_SPIN_RPM: f64 = 60_000.0;

/// Source-bounded operating inputs for the apparatus in US 1,102,653.
///
/// The facsimile supplies no absolute dimensions, masses, inertias, burn
/// rate, thrust, or spin speed. The two RPM values are therefore declared
/// visitor inputs. The kernel uses them only for torque-free rigid-body pose
/// and ideal instrument-isolation kinematics; it does not infer propulsion
/// performance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GoddardApparatusParams {
    /// Time at which the two torque-free display poses are sampled, in seconds.
    pub elapsed_seconds: f64,
    /// Declared primary-rocket spin speed, in revolutions per minute.
    pub primary_spin_rpm: f64,
    /// Declared gyroscope rotor speed, in revolutions per minute.
    pub gyro_spin_rpm: f64,
    /// Tapered-tube length divided by its longest diameter.
    pub tube_length_ratio: f64,
    /// Source-sequence animation coordinate: nested at 0, released at 1.
    pub auxiliary_release_fraction: f64,
    /// Whether the primary charge is in the source's "substantially consumed" state.
    pub primary_charge_substantially_consumed: bool,
    /// Whether the Claim 7 gyroscope restraint is present.
    pub gyro_enabled: bool,
}

/// Typed refusal from the source-bounded Goddard apparatus kernel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GoddardApparatusError {
    /// A named scalar input was NaN or infinite.
    NonFiniteInput(&'static str),
    /// A named scalar input was outside the bounded teaching domain.
    InputOutsideDomain(&'static str),
    /// The generic rigid-body owner refused the normalized pose calculation.
    RigidBody(DynamicsError),
}

impl From<DynamicsError> for GoddardApparatusError {
    fn from(value: DynamicsError) -> Self {
        Self::RigidBody(value)
    }
}

/// One source-bounded apparatus state derived from the generic rigid-body owner.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GoddardApparatusResult {
    /// Canonical `(w, x, y, z)` body-to-world pose of the primary rocket.
    pub primary_quaternion: [f64; 4],
    /// Canonical `(w, x, y, z)` body-to-world pose of the gyroscope rotor.
    pub gyro_quaternion: [f64; 4],
    /// Primary angular velocity in radians per second.
    pub primary_angular_velocity_rad_per_sec: f64,
    /// Gyroscope angular velocity in radians per second.
    pub gyro_angular_velocity_rad_per_sec: f64,
    /// Idealized instrument-support angular velocity in the world frame.
    pub camera_support_angular_velocity_rad_per_sec: f64,
    /// Primary rim speed per metre of radius, in metres per second per metre.
    pub primary_rim_speed_per_radius_mps_per_m: f64,
    /// The admitted tube length divided by longest diameter.
    pub tube_length_ratio: f64,
    /// Signed margin above or below Claim 2's `L / D = 3` lower bound.
    pub claim_2_ratio_margin: f64,
    /// Whether the Claim 2 tapered-tube ratio is satisfied.
    pub claim_2_satisfied: bool,
    /// Whether any requested auxiliary release follows substantial consumption.
    pub claim_1_sequence_satisfied: bool,
    /// Whether the auxiliary rocket is still nested in firing tube 24.
    pub auxiliary_nested: bool,
    /// Whether the Claim 7 gyroscope restraint is present.
    pub gyro_enabled: bool,
}

fn admit_apparatus_params(params: GoddardApparatusParams) -> Result<(), GoddardApparatusError> {
    for (name, value) in [
        ("elapsed_seconds", params.elapsed_seconds),
        ("primary_spin_rpm", params.primary_spin_rpm),
        ("gyro_spin_rpm", params.gyro_spin_rpm),
        ("tube_length_ratio", params.tube_length_ratio),
        (
            "auxiliary_release_fraction",
            params.auxiliary_release_fraction,
        ),
    ] {
        if !value.is_finite() {
            return Err(GoddardApparatusError::NonFiniteInput(name));
        }
    }

    if !(0.0..=GODDARD_MAX_ELAPSED_SECONDS).contains(&params.elapsed_seconds) {
        return Err(GoddardApparatusError::InputOutsideDomain("elapsed_seconds"));
    }
    if !(0.0..=GODDARD_MAX_PRIMARY_SPIN_RPM).contains(&params.primary_spin_rpm) {
        return Err(GoddardApparatusError::InputOutsideDomain(
            "primary_spin_rpm",
        ));
    }
    if !(0.0..=GODDARD_MAX_GYRO_SPIN_RPM).contains(&params.gyro_spin_rpm) {
        return Err(GoddardApparatusError::InputOutsideDomain("gyro_spin_rpm"));
    }
    if !(GODDARD_MIN_TUBE_LENGTH_RATIO..=GODDARD_MAX_TUBE_LENGTH_RATIO)
        .contains(&params.tube_length_ratio)
    {
        return Err(GoddardApparatusError::InputOutsideDomain(
            "tube_length_ratio",
        ));
    }
    if !(0.0..=1.0).contains(&params.auxiliary_release_fraction) {
        return Err(GoddardApparatusError::InputOutsideDomain(
            "auxiliary_release_fraction",
        ));
    }
    Ok(())
}

fn normalized_torque_free_pose(
    angular_velocity_body: Vec3,
    elapsed_seconds: f64,
) -> Result<UnitQuaternion, GoddardApparatusError> {
    if elapsed_seconds == 0.0 {
        return Ok(UnitQuaternion::IDENTITY);
    }

    // The patent gives no mass properties. A unit spherical inertia is used
    // only because, for torque-free spin, its orientation depends on the
    // declared angular velocity and elapsed time but not on an invented mass
    // or moment. No energy or force from this normalized body is published.
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0))?;
    let state = RigidBodyState::new(
        Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY)?,
        Vec3::ZERO,
        angular_velocity_body,
    )?;
    let receipt = RigidBodyIntegrator::new(Gravity::ZERO).step(
        state,
        properties,
        Wrench::ZERO,
        elapsed_seconds,
    )?;
    Ok(receipt.state_after.pose().orientation())
}

/// Steps the source-bounded US 1,102,653 teaching model.
///
/// `fs-mbd` owns both torque-free display poses. The remaining outputs are
/// source-sequence predicates or exact unit conversions. In particular this
/// function returns no liquid-propellant, de Laval, Mach, thrust, or trajectory
/// quantity because the patent supplies none of the inputs needed for them.
pub fn step_goddard_apparatus(
    params: &GoddardApparatusParams,
) -> Result<GoddardApparatusResult, GoddardApparatusError> {
    admit_apparatus_params(*params)?;

    let rpm_to_rad_per_sec = core::f64::consts::TAU / 60.0;
    let primary_angular_velocity_rad_per_sec = params.primary_spin_rpm * rpm_to_rad_per_sec;
    let gyro_angular_velocity_rad_per_sec = params.gyro_spin_rpm * rpm_to_rad_per_sec;
    let primary_quaternion = normalized_torque_free_pose(
        Vec3::new(0.0, primary_angular_velocity_rad_per_sec, 0.0),
        params.elapsed_seconds,
    )?
    .components();
    let gyro_quaternion = normalized_torque_free_pose(
        Vec3::new(gyro_angular_velocity_rad_per_sec, 0.0, 0.0),
        params.elapsed_seconds,
    )?
    .components();

    let auxiliary_nested = params.auxiliary_release_fraction == 0.0;
    Ok(GoddardApparatusResult {
        primary_quaternion,
        gyro_quaternion,
        primary_angular_velocity_rad_per_sec,
        gyro_angular_velocity_rad_per_sec,
        camera_support_angular_velocity_rad_per_sec: if params.gyro_enabled
            && params.gyro_spin_rpm > 0.0
        {
            0.0
        } else {
            primary_angular_velocity_rad_per_sec
        },
        primary_rim_speed_per_radius_mps_per_m: primary_angular_velocity_rad_per_sec,
        tube_length_ratio: params.tube_length_ratio,
        claim_2_ratio_margin: params.tube_length_ratio - GODDARD_CLAIM_2_MIN_TUBE_LENGTH_RATIO,
        claim_2_satisfied: params.tube_length_ratio >= GODDARD_CLAIM_2_MIN_TUBE_LENGTH_RATIO,
        claim_1_sequence_satisfied: auxiliary_nested
            || params.primary_charge_substantially_consumed,
        auxiliary_nested,
        gyro_enabled: params.gyro_enabled,
    })
}

#[derive(Debug, Clone)]
/// Inputs for the adjacent liquid-propellant rocket teaching model.
pub struct GoddardParams {
    /// Combustion-chamber pressure in pounds per square inch.
    pub chamber_pressure_psi: f64,
    /// Declared propellant mass flow in kilograms per second.
    pub fuel_flow_kg_per_sec: f64,
    /// Declared nozzle throat area in square centimetres.
    pub throat_area_cm2: f64,
    /// Declared exit-area to throat-area ratio.
    pub expansion_ratio: f64,
}

#[derive(Debug, Clone)]
/// Outputs from the adjacent liquid-propellant rocket teaching model.
pub struct GoddardResult {
    /// Echoed chamber pressure in pounds per square inch.
    pub chamber_pressure_psi: f64,
    /// Chamber pressure converted to pascals.
    pub chamber_pressure_pa: f64,
    /// Estimated ideal exhaust speed in metres per second.
    pub exhaust_velocity_mps: f64,
    /// Momentum-thrust estimate in newtons.
    pub thrust_newtons: f64,
    /// Estimated specific impulse in seconds.
    pub specific_impulse_sec: f64,
    /// Estimated nozzle-exit Mach number.
    pub mach_exit: f64,
}

/// Steps the adjacent liquid-propellant rocket teaching model.
///
/// This function is not the source owner for US 1,102,653; that record calls
/// [`step_goddard_apparatus`] instead.
pub fn step_goddard_rocket(params: &GoddardParams) -> GoddardResult {
    let chamber_pressure_pa = params.chamber_pressure_psi * 6894.76;
    let gamma = 1.24; // Combustion products heat capacity ratio
    let chamber_temp_k = (2400.0 + (chamber_pressure_pa / 2.4e6) * 400.0).round();
    let gas_constant_r = 365.0; // J/(kg*K) for gasoline + liquid O2
    let expansion = params.expansion_ratio.max(1.4);

    // Supersonic Mach number at exit via area-Mach relation
    let mach_exit =
        ((2.0 / (gamma - 1.0)) * (params.expansion_ratio.powf(2.0 / (gamma + 1.0)) - 1.0)).sqrt();
    let exhaust_velocity_mps = (((2.0 * gamma) / (gamma - 1.0))
        * gas_constant_r
        * chamber_temp_k
        * (1.0 - 1.0 / expansion.powf(gamma - 1.0)))
    .sqrt()
    .round();
    let thrust_newtons = (params.fuel_flow_kg_per_sec * exhaust_velocity_mps).round();
    let specific_impulse_sec = exhaust_velocity_mps / 9.80665;

    GoddardResult {
        chamber_pressure_psi: params.chamber_pressure_psi,
        chamber_pressure_pa,
        exhaust_velocity_mps,
        thrust_newtons,
        specific_impulse_sec: (specific_impulse_sec * 10.0).round() / 10.0,
        mach_exit: (mach_exit * 100.0).round() / 100.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn source_apparatus() -> GoddardApparatusParams {
        GoddardApparatusParams {
            elapsed_seconds: 0.25,
            primary_spin_rpm: 120.0,
            gyro_spin_rpm: 6_000.0,
            tube_length_ratio: 4.5,
            auxiliary_release_fraction: 0.0,
            primary_charge_substantially_consumed: false,
            gyro_enabled: true,
        }
    }

    #[test]
    fn apparatus_uses_rigid_body_poses_and_source_claim_predicates() {
        let result = step_goddard_apparatus(&source_apparatus()).expect("admitted apparatus");
        assert!(
            result
                .primary_quaternion
                .iter()
                .all(|value| value.is_finite())
        );
        assert!(result.gyro_quaternion.iter().all(|value| value.is_finite()));
        assert!(result.claim_2_satisfied);
        assert!(result.claim_1_sequence_satisfied);
        assert!(result.auxiliary_nested);
        assert_eq!(result.camera_support_angular_velocity_rad_per_sec, 0.0);
    }

    #[test]
    fn claim_two_break_probe_reports_the_signed_ratio_margin() {
        let mut params = source_apparatus();
        params.tube_length_ratio = 2.5;
        let result = step_goddard_apparatus(&params).expect("bounded break probe");
        assert!(!result.claim_2_satisfied);
        assert_eq!(result.claim_2_ratio_margin, -0.5);
    }

    #[test]
    fn auxiliary_sequence_and_gyro_omission_fail_visibly() {
        let mut params = source_apparatus();
        params.auxiliary_release_fraction = 0.5;
        params.gyro_enabled = false;
        let result = step_goddard_apparatus(&params).expect("bounded break probe");
        assert!(!result.claim_1_sequence_satisfied);
        assert!(!result.auxiliary_nested);
        assert_eq!(
            result.camera_support_angular_velocity_rad_per_sec,
            result.primary_angular_velocity_rad_per_sec
        );

        params.primary_charge_substantially_consumed = true;
        assert!(
            step_goddard_apparatus(&params)
                .expect("source sequence")
                .claim_1_sequence_satisfied
        );

        params.gyro_enabled = true;
        params.gyro_spin_rpm = 0.0;
        let stopped_gyro = step_goddard_apparatus(&params).expect("stopped gyro state");
        assert_eq!(
            stopped_gyro.camera_support_angular_velocity_rad_per_sec,
            stopped_gyro.primary_angular_velocity_rad_per_sec,
            "a present but stopped gyroscope cannot isolate the support",
        );
    }

    #[test]
    fn apparatus_refuses_non_finite_and_unbounded_inputs() {
        let mut params = source_apparatus();
        params.elapsed_seconds = f64::NAN;
        assert_eq!(
            step_goddard_apparatus(&params),
            Err(GoddardApparatusError::NonFiniteInput("elapsed_seconds"))
        );
        params.elapsed_seconds = 0.0;
        params.tube_length_ratio = GODDARD_MAX_TUBE_LENGTH_RATIO + 0.1;
        assert_eq!(
            step_goddard_apparatus(&params),
            Err(GoddardApparatusError::InputOutsideDomain(
                "tube_length_ratio"
            ))
        );
    }
}
