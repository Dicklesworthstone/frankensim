//! Prescribed-phase unbalance loads, including tangential angular acceleration.
//!
//! An eccentric material point has acceleration
//! `e*(-Omega²*e_r + alpha*e_theta)`. Its reaction on the nominal rotor is
//! therefore `U*(Omega²*e_r - alpha*e_theta)`, where U=m*e. Station moments
//! use the original transpose of `[x+z*theta_y, y-z*theta_x]`.
//! See Ansys Theory Reference, §14.4.2.4, Mass Unbalance Transient Forces:
//! <https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/ans_thry/thy_tool19.html>.
//!
//! This is the existing small-unbalance perturbation model. It does not add
//! another mass/inertia, supply a spin drive, or change the fixed-speed rotor
//! operator into a complete variable-spin or nonlinear bearing model.

use super::{
    MAX_ROTOR_STATIONS, RigidRotor, RotorError, RotorUnbalance, add_compensated, finite, poll,
};

/// Two complex physical moments of all original unbalance stations.
///
/// Preparation is O(stations); each load or pullback uses fixed-size work and
/// no allocation. Phase offsets and axial stations remain in the original
/// stationary rotor coordinates. No spectral or modal reduction is involved.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedRotorUnbalance {
    // [sum U cos(phi), sum U sin(phi), sum z U cos(phi), sum z U sin(phi)].
    moments: [f64; 4],
}

struct LoadSample {
    radial: [f64; 4],
    tangential: [f64; 4],
    force: [f64; 4],
}

/// Pullback of a scalar objective through the prescribed unbalance load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorUnbalanceKinematicsCotangent {
    /// Derivative with respect to the common shaft phase [objective/rad].
    pub phase: f64,
    /// Derivative with respect to signed shaft speed [objective/(rad/s)].
    pub spin: f64,
    /// Derivative with respect to signed angular acceleration [objective/(rad/s²)].
    pub angular_acceleration: f64,
}

impl RigidRotor {
    /// Prepare the original station geometry for transient unbalance forcing.
    ///
    /// Each stored phase is an offset from the common phase subsequently
    /// supplied to `force_at`. Unlike a steady harmonic force, that query
    /// includes the tangential reaction during startup, braking and reversal.
    /// Source moments must be finite and representable before speed scaling.
    /// Cancellation or invalid input returns no partial prepared source.
    pub fn prepare_unbalance(
        &self,
        unbalances: &[RotorUnbalance],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRotorUnbalance, RotorError> {
        poll(&mut cancelled)?;
        if unbalances.len() > MAX_ROTOR_STATIONS {
            return Err(RotorError::Invalid("rotor unbalance station count"));
        }
        let mut moments = [0.0; 4];
        let mut correction = [0.0; 4];
        for point in unbalances {
            poll(&mut cancelled)?;
            finite(point.station_m)?;
            finite(point.phase_rad)?;
            let u = point.mass_eccentricity_kg_m;
            if !u.is_finite() || u < 0.0 {
                return Err(RotorError::Invalid(
                    "nonnegative finite rotor mass-eccentricity required",
                ));
            }
            let x = product(u, fs_math::det::cos(point.phase_rad))?;
            let y = product(u, fs_math::det::sin(point.phase_rad))?;
            let terms = [
                x,
                y,
                product(point.station_m, x)?,
                product(point.station_m, y)?,
            ];
            for i in 0..4 {
                add_compensated(&mut moments[i], &mut correction[i], terms[i])?;
            }
        }
        for i in 0..4 {
            moments[i] = finite(moments[i] + correction[i])?;
        }
        poll(&mut cancelled)?;
        Ok(PreparedRotorUnbalance { moments })
    }
}

impl PreparedRotorUnbalance {
    /// Original generalized force `[Fx,Fy,tau_x,tau_y]` in `[N,N,Nm,Nm]`.
    ///
    /// All stations share the supplied continuous shaft phase [rad], signed
    /// speed [rad/s] and angular acceleration [rad/s²]. Phase is independent
    /// of speed here: callers integrate their physical speed history rather
    /// than forming `current_speed * absolute_time`. A timestep owner should
    /// sample these quantities at its actual force-evaluation time.
    ///
    /// The source evaluates `Omega²*radial - alpha*tangential`. In particular,
    /// zero speed with nonzero angular acceleration generally has nonzero
    /// force. The retained nominal rotor inertia is unchanged. This is a load
    /// source, not a spin trajectory, shaft-drive energy audit or state-error
    /// certificate. Finite cancellation checks surround the bounded query.
    pub fn force_at(
        &self,
        phase_rad: f64,
        spin_rad_s: f64,
        angular_acceleration_rad_s2: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 4], RotorError> {
        poll(&mut cancelled)?;
        let sample = self.evaluate(phase_rad, spin_rad_s, angular_acceleration_rad_s2)?;
        poll(&mut cancelled)?;
        Ok(sample.force)
    }

    /// Exact derivative of the admitted physical load formula at fixed source
    /// geometry. The seed is conjugate to `[N,N,Nm,Nm]`. Compose its three
    /// outputs with the caller's phase/speed/acceleration history derivatives.
    /// This does not differentiate rounding, admission, or a rotor solver.
    pub fn force_vjp(
        &self,
        phase_rad: f64,
        spin_rad_s: f64,
        angular_acceleration_rad_s2: f64,
        seed: [f64; 4],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorUnbalanceKinematicsCotangent, RotorError> {
        poll(&mut cancelled)?;
        let LoadSample {
            radial, tangential, ..
        } = self.evaluate(phase_rad, spin_rad_s, angular_acceleration_rad_s2)?;
        let mut sums = [0.0; 3];
        let mut correction = [0.0; 3];
        for i in 0..4 {
            finite(seed[i])?;
            let phase = finite(
                product(product(tangential[i], spin_rad_s)?, spin_rad_s)?
                    + product(radial[i], angular_acceleration_rad_s2)?,
            )?;
            let spin = product(2.0, product(radial[i], spin_rad_s)?)?;
            for (j, derivative) in [phase, spin, -tangential[i]].into_iter().enumerate() {
                add_compensated(
                    &mut sums[j],
                    &mut correction[j],
                    product(seed[i], derivative)?,
                )?;
            }
        }
        for i in 0..3 {
            sums[i] = finite(sums[i] + correction[i])?;
        }
        poll(&mut cancelled)?;
        Ok(RotorUnbalanceKinematicsCotangent {
            phase: sums[0],
            spin: sums[1],
            angular_acceleration: sums[2],
        })
    }

    fn evaluate(
        &self,
        phase_rad: f64,
        spin_rad_s: f64,
        angular_acceleration_rad_s2: f64,
    ) -> Result<LoadSample, RotorError> {
        finite(spin_rad_s)?;
        finite(angular_acceleration_rad_s2)?;
        let (radial, tangential) = self.directions(phase_rad)?;
        let mut force = [0.0; 4];
        for i in 0..4 {
            // Do not form Omega² first: it can overflow even when U*Omega²
            // is representable. This order has monotone magnitude scaling.
            let centrifugal = product(product(radial[i], spin_rad_s)?, spin_rad_s)?;
            force[i] = finite(centrifugal - product(tangential[i], angular_acceleration_rad_s2)?)?;
        }
        Ok(LoadSample {
            radial,
            tangential,
            force,
        })
    }

    fn directions(&self, phase_rad: f64) -> Result<([f64; 4], [f64; 4]), RotorError> {
        finite(phase_rad)?;
        let (sin, cos) = (fs_math::det::sin(phase_rad), fs_math::det::cos(phase_rad));
        let mut rotated = [0.0; 4];
        for i in [0, 2] {
            // Admit the phase-projected moments before speed scaling as well;
            // a vanished product could become a resolved load at high speed.
            product(cos, self.moments[i])?;
            let y_sin = product(sin, self.moments[i + 1])?;
            product(sin, self.moments[i])?;
            let y_cos = product(cos, self.moments[i + 1])?;
            rotated[i] = finite(cos.mul_add(self.moments[i], -y_sin))?;
            rotated[i + 1] = finite(sin.mul_add(self.moments[i], y_cos))?;
        }
        let [x, y, zx, zy] = rotated;
        Ok(([x, y, -zy, zx], [-y, x, -zx, -zy]))
    }
}

fn product(a: f64, b: f64) -> Result<f64, RotorError> {
    let value = finite(a * b)?;
    if a != 0.0 && b != 0.0 && value == 0.0 {
        Err(RotorError::Invalid(
            "unrepresentable rotor unbalance product",
        ))
    } else {
        Ok(value)
    }
}
