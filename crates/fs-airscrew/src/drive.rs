//! A declared engine rating and passive transmission on one propeller coordinate.
//!
//! This is the existing flat-torque/constant-power approximation, made explicit
//! and selectable. The source supplies mechanical shaft power; fuel chemistry,
//! electrical conversion, measured engine maps and thermal derating are not
//! inferred. The common output coordinate may drive several gear-locked props;
//! its torque is the TOTAL drive torque, not a per-propeller allocation.
//! Both engine and propeller rates are relative to the carrier/crankcase;
//! aerodynamic shaft work uses absolute spin and remains a separate port.

use super::{Refusal, refuse};
use fs_blake3::hash_domain;

/// Prepared mechanical source and transmission, in SI units.
///
/// The engine turns `engine_to_prop_ratio` times as fast as the common prop
/// coordinate. Before losses, torque is constant below the rated speed and
/// power is constant above it. A command scales both source limits. Efficiency
/// applies only to transmitted source power; bearing drag remains when the
/// engine is off. Parameters are declared estimates unless independently sourced.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShaftDrive {
    rated_power_w: f64,
    rated_engine_speed_rad_s: f64,
    engine_to_prop_ratio: f64,
    transmission_efficiency: f64,
    prop_viscous_drag_nms: f64,
    rated_prop_speed_rad_s: f64,
    prop_torque_limit_nm: f64,
}

/// The instantaneous source, delivered drive and physical loss channels.
///
/// All powers use the same supplied shaft speed. These are source/transmission
/// powers, not propeller aerodynamic power or a whole-aircraft energy balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DriveOperatingPoint {
    /// Nonnegative engine speed relative to the crankcase [rad/s].
    pub engine_speed_rad_s: f64,
    /// Mechanical torque supplied by the engine before transmission loss [N m].
    pub engine_torque_nm: f64,
    /// Engine torque times engine speed [W], before transmission losses.
    pub engine_power_w: f64,
    /// Total net drive torque conjugate to the common propeller angle [N m].
    /// Negative values brake the shaft; they are never clipped to zero.
    pub prop_torque_nm: f64,
    /// Net drive torque times common propeller speed [W].
    pub prop_power_w: f64,
    /// Nonnegative loss of the declared constant-efficiency transmission [W].
    pub transmission_heat_w: f64,
    /// Nonnegative viscous bearing loss, `b * Omega_prop^2` [W].
    pub bearing_heat_w: f64,
    /// Engine minus net propeller power and the two physical losses [W].
    /// Floating-point closure error stays separate from heat.
    pub power_residual_w: f64,
}

impl ShaftDrive {
    /// Prepare a positive rating and speed ratio, efficiency in `[0,1]`, and
    /// nonnegative viscous drag on the TOTAL common propeller coordinate.
    ///
    /// # Errors
    /// `shaft-drive-invalid` for nonfinite or unrepresentable parameters.
    pub fn new(
        rated_power_w: f64,
        rated_engine_speed_rad_s: f64,
        engine_to_prop_ratio: f64,
        transmission_efficiency: f64,
        prop_viscous_drag_nms: f64,
    ) -> Result<Self, Refusal> {
        if ![
            rated_power_w,
            rated_engine_speed_rad_s,
            engine_to_prop_ratio,
            transmission_efficiency,
            prop_viscous_drag_nms,
        ]
        .into_iter()
        .all(f64::is_finite)
            || rated_power_w <= 0.0
            || rated_engine_speed_rad_s <= 0.0
            || engine_to_prop_ratio <= 0.0
            || !(0.0..=1.0).contains(&transmission_efficiency)
            || prop_viscous_drag_nms < 0.0
        {
            return Err(invalid(
                "invalid source rating, gear ratio or passive loss law",
            ));
        }
        let rated_prop_speed_rad_s = rated_engine_speed_rad_s * engine_to_prop_ratio.recip();
        let prop_torque_limit_nm = rated_power_w / rated_prop_speed_rad_s;
        let engine_torque_limit_nm = rated_power_w / rated_engine_speed_rad_s;
        if !rated_prop_speed_rad_s.is_finite()
            || rated_prop_speed_rad_s <= 0.0
            || !prop_torque_limit_nm.is_finite()
            || prop_torque_limit_nm <= 0.0
            || !engine_torque_limit_nm.is_finite()
            || engine_torque_limit_nm <= 0.0
        {
            return Err(invalid("unrepresentable reflected shaft rating"));
        }
        Ok(Self {
            rated_power_w,
            rated_engine_speed_rad_s,
            engine_to_prop_ratio,
            transmission_efficiency: if transmission_efficiency == 0.0 {
                0.0
            } else {
                transmission_efficiency
            },
            prop_viscous_drag_nms: if prop_viscous_drag_nms == 0.0 {
                0.0
            } else {
                prop_viscous_drag_nms
            },
            rated_prop_speed_rad_s,
            prop_torque_limit_nm,
        })
    }

    /// Original 12 hp / 1025 rpm / 23:8 chain approximation, with the original
    /// ideal transmission. No new historical loss or control setting is assumed.
    #[must_use]
    pub fn wright_reference() -> Self {
        Self::new(
            8948.0,
            1025.0 / 60.0 * 2.0 * core::f64::consts::PI,
            23.0 / 8.0,
            1.0,
            0.0,
        )
        .expect("finite positive reference drive")
    }

    /// Source rating before transmission loss [W].
    #[must_use]
    pub const fn rated_power_w(&self) -> f64 {
        self.rated_power_w
    }

    /// Common-shaft speed at the torque/power transition [rad/s].
    #[must_use]
    pub const fn rated_prop_speed_rad_s(&self) -> f64 {
        self.rated_prop_speed_rad_s
    }

    pub(super) fn unscaled_prop_torque(&self, omega_prop_rad_s: f64) -> f64 {
        if omega_prop_rad_s <= self.rated_prop_speed_rad_s {
            self.prop_torque_limit_nm
        } else {
            self.rated_power_w / omega_prop_rad_s
        }
    }

    /// Evaluate a nonnegative common-shaft rate relative to the carrier and
    /// a source command `[0,1]`.
    ///
    /// Zero command removes engine torque, not bearing resistance. Zero speed
    /// has finite starting torque and zero power. The law is unidirectional:
    /// reverse rotation refuses instead of inventing an engine or generator map.
    /// Constant efficiency is a declared approximation and does not represent
    /// clutch slip, temperature feedback, static friction or an engine redline.
    /// Time integrators must evaluate this law at their actual stage states;
    /// holding one sampled torque does not enforce a nonlinear power cap over
    /// an arbitrary finite acceleration interval.
    ///
    /// # Errors
    /// `shaft-drive-invalid` for command/speed domain or arithmetic overflow.
    pub fn evaluate(
        &self,
        omega_prop_rad_s: f64,
        command: f64,
    ) -> Result<DriveOperatingPoint, Refusal> {
        if !omega_prop_rad_s.is_finite()
            || omega_prop_rad_s < 0.0
            || !command.is_finite()
            || !(0.0..=1.0).contains(&command)
        {
            return Err(invalid(
                "require finite nonnegative shaft speed and command in [0,1]",
            ));
        }
        let source_prop_torque = command * self.unscaled_prop_torque(omega_prop_rad_s);
        let engine_speed = omega_prop_rad_s * self.engine_to_prop_ratio;
        let engine_torque = source_prop_torque / self.engine_to_prop_ratio;
        let engine_power = engine_torque * engine_speed;
        let bearing_torque = self.prop_viscous_drag_nms * omega_prop_rad_s;
        let prop_torque = self.transmission_efficiency * source_prop_torque - bearing_torque;
        let prop_power = prop_torque * omega_prop_rad_s;
        let transmission_heat = (1.0 - self.transmission_efficiency) * engine_power;
        let bearing_heat = bearing_torque * omega_prop_rad_s;
        let residual = ((engine_power - prop_power) - transmission_heat) - bearing_heat;
        if ![
            engine_speed,
            engine_torque,
            engine_power,
            prop_torque,
            prop_power,
            transmission_heat,
            bearing_heat,
            residual,
        ]
        .into_iter()
        .all(f64::is_finite)
            || (omega_prop_rad_s > 0.0 && engine_speed == 0.0)
            || (command > 0.0
                && (source_prop_torque == 0.0
                    || engine_torque == 0.0
                    || (omega_prop_rad_s > 0.0 && engine_power == 0.0)))
        {
            return Err(invalid(
                "unrepresentable shaft speed, torque or source power",
            ));
        }
        Ok(DriveOperatingPoint {
            engine_speed_rad_s: engine_speed,
            engine_torque_nm: engine_torque,
            engine_power_w: engine_power,
            prop_torque_nm: prop_torque,
            prop_power_w: prop_power,
            transmission_heat_w: transmission_heat,
            bearing_heat_w: bearing_heat,
            power_residual_w: residual,
        })
    }

    /// Identity of the declared physical parameters, independent of preparation.
    #[must_use]
    pub fn digest(&self) -> String {
        let mut bytes = [0; 40];
        for (slot, value) in bytes.chunks_exact_mut(8).zip([
            self.rated_power_w,
            self.rated_engine_speed_rad_s,
            self.engine_to_prop_ratio,
            self.transmission_efficiency,
            self.prop_viscous_drag_nms,
        ]) {
            slot.copy_from_slice(&value.to_bits().to_le_bytes());
        }
        hash_domain("org.frankensim.fs-airscrew.shaft-drive.v1", &bytes).to_hex()
    }
}

fn invalid(message: &str) -> Refusal {
    refuse(
        "shaft-drive-invalid",
        message.into(),
        "use a finite positive drive rating/gear ratio, passive losses and a command in [0,1]",
    )
}
