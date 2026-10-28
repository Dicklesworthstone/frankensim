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

/// Exact conventional pound-force per square inch in SI pressure units.
pub const GODDARD_PSI_TO_PA: f64 = 6_894.757_293_168_361;

/// Largest exit/throat area ratio admitted by this bounded teaching model.
pub const GODDARD_MAX_EXPANSION_RATIO: f64 = 1_000.0;

/// Inputs for an ideal calorically-perfect rocket nozzle.
///
/// Temperature, gas constant and heat-capacity ratio are independent declared
/// thermodynamic assumptions, not inferred from chamber pressure or attributed
/// to a historical propellant. The throat determines mass flow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GoddardParams {
    /// Absolute chamber stagnation pressure in pounds per square inch.
    pub chamber_pressure_psi: f64,
    /// Chamber stagnation temperature in kelvin.
    pub chamber_temperature_k: f64,
    /// Frozen ideal-gas specific gas constant in J/(kg K).
    pub specific_gas_constant_j_kg_k: f64,
    /// Constant cp/cv, greater than one and no greater than 5/3.
    pub heat_capacity_ratio: f64,
    /// Declared nozzle throat area in square centimetres.
    pub throat_area_cm2: f64,
    /// Exit/throat area ratio in [1, 1000]; one is a sonic converging exit.
    pub expansion_ratio: f64,
    /// Absolute external pressure in pascals; zero explicitly means vacuum.
    pub ambient_pressure_pa: f64,
}

/// Ideal-nozzle outputs in SI, without presentation rounding.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GoddardResult {
    /// Echoed chamber pressure in pounds per square inch.
    pub chamber_pressure_psi: f64,
    /// Chamber pressure converted to pascals.
    pub chamber_pressure_pa: f64,
    /// Static nozzle-exit temperature in kelvin.
    pub exit_temperature_k: f64,
    /// Static nozzle-exit pressure in pascals.
    pub exit_pressure_pa: f64,
    /// Choked throat mass flow in kilograms per second.
    pub mass_flow_kg_per_sec: f64,
    /// Ideal nozzle-exit gas speed in metres per second.
    pub exhaust_velocity_mps: f64,
    /// Exhaust momentum flux in newtons.
    pub momentum_thrust_newtons: f64,
    /// Exit pressure minus ambient pressure, times exit area, in newtons.
    pub pressure_thrust_newtons: f64,
    /// Momentum plus pressure thrust in newtons.
    pub thrust_newtons: f64,
    /// Total thrust divided by mass flow and standard gravity, in seconds.
    pub specific_impulse_sec: f64,
    /// Chamber pressure times throat area divided by mass flow, in m/s.
    pub characteristic_velocity_mps: f64,
    /// Supersonic nozzle-exit Mach number, or exactly one for unit area ratio.
    pub mach_exit: f64,
}

/// Explicit refusal rather than a fabricated gas state or a clamped geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GoddardRocketError {
    /// A named scalar input is not finite.
    NonFiniteInput(&'static str),
    /// A named input is outside this model's physical or work bounds.
    InputOutsideDomain(&'static str),
    /// Overexpanded operation needs a shock/separation model not supplied here.
    UnsupportedBackPressure {
        /// Requested ambient pressure [Pa].
        ambient_pressure_pa: f64,
        /// Isentropic exit pressure [Pa].
        exit_pressure_pa: f64,
    },
    /// A necessary positive finite quantity cannot be represented.
    UnrepresentableOutput(&'static str),
    /// The bounded supersonic area-Mach inversion did not converge.
    NoConvergence,
}

impl core::fmt::Display for GoddardRocketError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ideal rocket nozzle: {self:?}")
    }
}
impl std::error::Error for GoddardRocketError {}

// Recover the rounding lost when forming 1+x. In particular gamma-1 may be
// small even though the resulting isentropic exponent is order one.
fn nozzle_log1p(x: f64) -> f64 {
    let y = 1.0 + x;
    fs_math::det::ln(y) + (x - (y - 1.0)) / y
}

fn nozzle_log1pmx_small(x: f64) -> f64 {
    let mut power = x * x;
    let mut sum = -0.5 * power;
    for n in 3..=24 {
        power *= -x;
        sum -= power / f64::from(n);
    }
    sum
}

fn nozzle_log_area(mach: f64, gamma: f64) -> f64 {
    let y = (mach - 1.0) * (mach + 1.0);
    let a = (gamma - 1.0) / (gamma + 1.0);
    // ln(A/A*) = [ln(1+a*y)/a - ln(1+y)]/2. Near choking,
    // cancel the equal linear terms analytically, not by floating subtraction.
    if y < 0.25 {
        0.5 * (nozzle_log1pmx_small(a * y) / a - nozzle_log1pmx_small(y))
    } else {
        0.5 * (nozzle_log1p(a * y) / a - nozzle_log1p(y))
    }
}

fn nozzle_supersonic_mach(area_ratio: f64, gamma: f64) -> Result<f64, GoddardRocketError> {
    if area_ratio == 1.0 {
        return Ok(1.0);
    }
    let target = fs_math::det::ln(area_ratio);
    // M=64 brackets A/A* <= 1000 for every admitted 1 < gamma <= 5/3.
    // Fixed bisection is allocation-free and gives a bounded work cost.
    let (mut lower, mut upper) = (1.0, 64.0);
    for _ in 0..80 {
        let mach = lower + 0.5 * (upper - lower);
        let residual = nozzle_log_area(mach, gamma) - target;
        if upper - lower <= 8.0 * f64::EPSILON * mach {
            return if residual.abs() <= 1e-12 {
                Ok(mach)
            } else {
                Err(GoddardRocketError::NoConvergence)
            };
        }
        if residual > 0.0 {
            upper = mach;
        } else {
            lower = mach;
        }
    }
    Err(GoddardRocketError::NoConvergence)
}

/// Resolve a quasi-steady, ideal calorically-perfect rocket nozzle.
///
/// This function is not the source owner for US 1,102,653; that record calls
/// [`step_goddard_apparatus`] instead. NASA Glenn's isentropic area-Mach and
/// rocket thrust relations are used without a propellant-specific fit:
/// <https://www.grc.nasa.gov/www/BGH/rktthsum.html>.
///
/// Mass flow is set by the sonic throat. The supersonic exit conserves that
/// flow and stagnation enthalpy; thrust is `mdot*ve + (pe-pa)*Ae`. This lane
/// admits pressure-matched and underexpanded operation, including vacuum.
/// Overexpanded operation refuses because internal shocks and separation
/// cannot be inferred from area ratio alone. No combustion chemistry,
/// variable heat capacities, nozzle loss, vehicle mass or trajectory is modeled.
pub fn step_goddard_rocket(params: &GoddardParams) -> Result<GoddardResult, GoddardRocketError> {
    for (name, value) in [
        ("chamber_pressure_psi", params.chamber_pressure_psi),
        ("chamber_temperature_k", params.chamber_temperature_k),
        (
            "specific_gas_constant_j_kg_k",
            params.specific_gas_constant_j_kg_k,
        ),
        ("heat_capacity_ratio", params.heat_capacity_ratio),
        ("throat_area_cm2", params.throat_area_cm2),
        ("expansion_ratio", params.expansion_ratio),
        ("ambient_pressure_pa", params.ambient_pressure_pa),
    ] {
        if !value.is_finite() {
            return Err(GoddardRocketError::NonFiniteInput(name));
        }
        let valid = match name {
            "heat_capacity_ratio" => value > 1.0 && value <= 5.0 / 3.0,
            "expansion_ratio" => (1.0..=GODDARD_MAX_EXPANSION_RATIO).contains(&value),
            "ambient_pressure_pa" => value >= 0.0,
            _ => value > 0.0,
        };
        if !valid {
            return Err(GoddardRocketError::InputOutsideDomain(name));
        }
    }
    let gamma = params.heat_capacity_ratio;
    let delta_gamma = gamma - 1.0;
    let chamber_pressure_pa = params.chamber_pressure_psi * GODDARD_PSI_TO_PA;
    let throat_area_m2 = params.throat_area_cm2 * 1e-4;
    let exit_area_m2 = throat_area_m2 * params.expansion_ratio;
    let rt0 = params.specific_gas_constant_j_kg_k * params.chamber_temperature_k;
    for (name, value) in [
        ("chamber_pressure_pa", chamber_pressure_pa),
        ("throat_area_m2", throat_area_m2),
        ("exit_area_m2", exit_area_m2),
        ("specific_gas_constant_times_temperature", rt0),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(GoddardRocketError::UnrepresentableOutput(name));
        }
    }
    let mach_exit = nozzle_supersonic_mach(params.expansion_ratio, gamma)?;
    let thermal_factor = 0.5 * delta_gamma * mach_exit * mach_exit;
    let exit_temperature_k = params.chamber_temperature_k / (1.0 + thermal_factor);
    let exit_pressure_pa = chamber_pressure_pa
        * fs_math::det::exp(-gamma / delta_gamma * nozzle_log1p(thermal_factor));
    if params.ambient_pressure_pa > exit_pressure_pa {
        return Err(GoddardRocketError::UnsupportedBackPressure {
            ambient_pressure_pa: params.ambient_pressure_pa,
            exit_pressure_pa,
        });
    }
    let sonic_factor =
        fs_math::det::exp(-0.5 * (gamma + 1.0) / delta_gamma * nozzle_log1p(0.5 * delta_gamma));
    let mass_flow_kg_per_sec = throat_area_m2
        * (chamber_pressure_pa / fs_math::det::sqrt(rt0))
        * fs_math::det::sqrt(gamma)
        * sonic_factor;
    let exhaust_velocity_mps = mach_exit
        * fs_math::det::sqrt(gamma * params.specific_gas_constant_j_kg_k * exit_temperature_k);
    let momentum_thrust_newtons = mass_flow_kg_per_sec * exhaust_velocity_mps;
    let pressure_thrust_newtons = (exit_pressure_pa - params.ambient_pressure_pa) * exit_area_m2;
    let thrust_newtons = momentum_thrust_newtons + pressure_thrust_newtons;
    let specific_impulse_sec = (thrust_newtons / mass_flow_kg_per_sec) / 9.80665;
    let characteristic_velocity_mps =
        fs_math::det::sqrt(rt0) / (fs_math::det::sqrt(gamma) * sonic_factor);
    for (name, value) in [
        ("exit_temperature_k", exit_temperature_k),
        ("exit_pressure_pa", exit_pressure_pa),
        ("mass_flow_kg_per_sec", mass_flow_kg_per_sec),
        ("exhaust_velocity_mps", exhaust_velocity_mps),
        ("momentum_thrust_newtons", momentum_thrust_newtons),
        ("thrust_newtons", thrust_newtons),
        ("specific_impulse_sec", specific_impulse_sec),
        ("characteristic_velocity_mps", characteristic_velocity_mps),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(GoddardRocketError::UnrepresentableOutput(name));
        }
    }
    if !pressure_thrust_newtons.is_finite() || pressure_thrust_newtons < 0.0 {
        return Err(GoddardRocketError::UnrepresentableOutput(
            "pressure_thrust_newtons",
        ));
    }
    Ok(GoddardResult {
        chamber_pressure_psi: params.chamber_pressure_psi,
        chamber_pressure_pa,
        exit_temperature_k,
        exit_pressure_pa,
        mass_flow_kg_per_sec,
        exhaust_velocity_mps,
        momentum_thrust_newtons,
        pressure_thrust_newtons,
        thrust_newtons,
        specific_impulse_sec,
        characteristic_velocity_mps,
        mach_exit,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ideal_nozzle() -> GoddardParams {
        GoddardParams {
            chamber_pressure_psi: 2.0e6 / GODDARD_PSI_TO_PA,
            chamber_temperature_k: 1800.0,
            specific_gas_constant_j_kg_k: 300.0,
            heat_capacity_ratio: 1.4,
            throat_area_cm2: 10.0,
            expansion_ratio: 27.0 / 16.0,
            ambient_pressure_pa: 0.0,
        }
    }

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 5e-12 * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn g1_mach_two_nozzle_closes_mass_energy_and_thrust_against_analytic_values() {
        let params = ideal_nozzle();
        let result = step_goddard_rocket(&params).unwrap();
        // For gamma=7/5 and M=2, A/A*=27/16 and T/T0=5/9.
        let exit_pressure = 2.0e6 * (5.0_f64 / 9.0).powf(3.5);
        let velocity = 2.0 * (1.4_f64 * 300.0 * 1000.0).sqrt();
        let rate = 2000.0 / (300.0_f64 * 1800.0).sqrt() * 1.4_f64.sqrt() * (5.0_f64 / 6.0).powi(3);
        near(result.mach_exit, 2.0);
        near(result.exit_temperature_k, 1000.0);
        near(result.exit_pressure_pa, exit_pressure);
        near(result.exhaust_velocity_mps, velocity);
        near(result.mass_flow_kg_per_sec, rate);
        let area = 0.001 * params.expansion_ratio;
        near(
            result.mass_flow_kg_per_sec,
            exit_pressure / 300000.0 * velocity * area,
        );
        near(
            0.5 * velocity * velocity,
            1050.0 * (1800.0 - result.exit_temperature_k),
        );
        near(
            result.thrust_newtons,
            rate * velocity + exit_pressure * area,
        );
        near(
            result.specific_impulse_sec,
            result.thrust_newtons / rate / 9.80665,
        );
        near(result.characteristic_velocity_mps, 2000.0 / rate);
        assert!(result.specific_impulse_sec > velocity / 9.80665);
        assert_eq!(step_goddard_rocket(&params).unwrap(), result);
    }

    #[test]
    fn sonic_limit_and_near_sonic_branch_do_not_use_a_false_mach_approximation() {
        let mut params = ideal_nozzle();
        params.expansion_ratio = 1.0;
        let sonic = step_goddard_rocket(&params).unwrap();
        assert_eq!(sonic.mach_exit, 1.0);
        near(sonic.exit_temperature_k, 1500.0);
        near(sonic.exit_pressure_pa, 2.0e6 * (5.0_f64 / 6.0).powf(3.5));
        params.expansion_ratio = f64::from_bits(1.0_f64.to_bits() + 1);
        let near_sonic = step_goddard_rocket(&params).unwrap();
        assert!(near_sonic.mach_exit > 1.0);
        near(
            near_sonic.mach_exit,
            1.0 + (1.2 * (params.expansion_ratio - 1.0)).sqrt(),
        );
    }

    #[test]
    fn throat_pressure_and_temperature_have_their_distinct_physical_scalings() {
        let p = ideal_nozzle();
        let base = step_goddard_rocket(&p).unwrap();
        for doubled in [
            GoddardParams {
                throat_area_cm2: 2.0 * p.throat_area_cm2,
                ..p
            },
            GoddardParams {
                chamber_pressure_psi: 2.0 * p.chamber_pressure_psi,
                ..p
            },
        ] {
            let result = step_goddard_rocket(&doubled).unwrap();
            near(result.mass_flow_kg_per_sec, 2.0 * base.mass_flow_kg_per_sec);
            near(result.thrust_newtons, 2.0 * base.thrust_newtons);
            near(result.exhaust_velocity_mps, base.exhaust_velocity_mps);
            near(result.specific_impulse_sec, base.specific_impulse_sec);
            assert_eq!(result.mach_exit, base.mach_exit);
        }
        let hot = step_goddard_rocket(&GoddardParams {
            chamber_temperature_k: 2.0 * p.chamber_temperature_k,
            ..p
        })
        .unwrap();
        near(
            hot.mass_flow_kg_per_sec,
            base.mass_flow_kg_per_sec / 2.0_f64.sqrt(),
        );
        near(
            hot.exhaust_velocity_mps,
            base.exhaust_velocity_mps * 2.0_f64.sqrt(),
        );
        near(hot.thrust_newtons, base.thrust_newtons);
    }

    #[test]
    fn ambient_pressure_changes_pressure_thrust_without_changing_choked_flow() {
        let p = ideal_nozzle();
        let vacuum = step_goddard_rocket(&p).unwrap();
        for ambient in [101325.0, vacuum.exit_pressure_pa] {
            let result = step_goddard_rocket(&GoddardParams {
                ambient_pressure_pa: ambient,
                ..p
            })
            .unwrap();
            assert_eq!(result.mass_flow_kg_per_sec, vacuum.mass_flow_kg_per_sec);
            assert_eq!(result.exhaust_velocity_mps, vacuum.exhaust_velocity_mps);
            near(
                vacuum.thrust_newtons - result.thrust_newtons,
                ambient * 0.001 * p.expansion_ratio,
            );
            if ambient == vacuum.exit_pressure_pa {
                assert_eq!(result.pressure_thrust_newtons, 0.0);
                near(
                    result.specific_impulse_sec,
                    result.exhaust_velocity_mps / 9.80665,
                );
            }
        }
        assert!(matches!(
            step_goddard_rocket(&GoddardParams {
                ambient_pressure_pa: vacuum.exit_pressure_pa + 1.0,
                ..p
            }),
            Err(GoddardRocketError::UnsupportedBackPressure { .. })
        ));
    }

    #[test]
    fn gamma_near_one_and_large_area_ratio_have_finite_physical_limits() {
        let p = GoddardParams {
            heat_capacity_ratio: f64::from_bits(1.0_f64.to_bits() + 1),
            expansion_ratio: 1000.0,
            ..ideal_nozzle()
        };
        let result = step_goddard_rocket(&p).unwrap();
        let m = result.mach_exit;
        // Independent gamma->1 limit: A/A*=exp((M^2-1)/2)/M.
        near(((m * m - 1.0) / 2.0).exp() / m, p.expansion_ratio);
        near(
            result.exit_pressure_pa,
            result.chamber_pressure_pa * (-m * m / 2.0).exp(),
        );
        let mono = step_goddard_rocket(&GoddardParams {
            heat_capacity_ratio: 5.0 / 3.0,
            ..p
        })
        .unwrap();
        assert!(mono.mach_exit > m && mono.mach_exit < 64.0);
    }

    #[test]
    fn invalid_thermodynamics_geometry_and_unrepresentable_states_refuse() {
        let p = ideal_nozzle();
        for invalid in [
            GoddardParams {
                chamber_pressure_psi: f64::NAN,
                ..p
            },
            GoddardParams {
                chamber_temperature_k: 0.0,
                ..p
            },
            GoddardParams {
                specific_gas_constant_j_kg_k: -1.0,
                ..p
            },
            GoddardParams {
                heat_capacity_ratio: 1.0,
                ..p
            },
            GoddardParams {
                heat_capacity_ratio: 1.7,
                ..p
            },
            GoddardParams {
                throat_area_cm2: 0.0,
                ..p
            },
            GoddardParams {
                expansion_ratio: 0.999,
                ..p
            },
            GoddardParams {
                expansion_ratio: 1001.0,
                ..p
            },
            GoddardParams {
                ambient_pressure_pa: -1.0,
                ..p
            },
        ] {
            assert!(step_goddard_rocket(&invalid).is_err());
        }
        assert!(matches!(
            step_goddard_rocket(&GoddardParams {
                throat_area_cm2: f64::from_bits(1),
                ..p
            }),
            Err(GoddardRocketError::UnrepresentableOutput("throat_area_m2"))
        ));
    }

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
