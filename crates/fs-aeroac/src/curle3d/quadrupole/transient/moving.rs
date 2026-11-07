//! Moving Lighthill stress moments on their original emission clock.
//!
//! The wave equation gives `p = partial_i partial_j[Q_ij(tau)/(4 pi R kappa)]`,
//! `kappa = 1 - n.v/c`. Both spatial derivatives act on the retarded clock AND
//! the convective Green amplitude. Replacing the stationary delay alone loses
//! these terms. Source strength remains the physical stress moment [N m]; no
//! density or extra Doppler multiplier belongs in it.
//! Tensor components and their time derivatives are expressed in fixed axes of
//! the medium rest frame. Uniform translation does not rotate that tensor;
//! rotating-body producers must supply the corresponding rotated derivatives.
//!
//! For uniform velocity, put `r = x-y(t)`, `M=v/c`, `b=1-|M|^2`, and
//! `D = sqrt(b r.r + (M.r)^2) = R kappa`. The derivatives of D and of
//! `tau=t-(M.r+D)/(c b)` give the three tensor coefficients below. This is the
//! moving-source form of the Lighthill Green solution, Lighthill (1952),
//! <https://doi.org/10.1098/rspa.1952.0060>. The retarded Jacobian follows
//! Farassat, NASA/TM-2007-214853,
//! <https://ntrs.nasa.gov/citations/20070010579>.
//!
//! Accelerating motion uses derivatives of the SAME retarded Green potential,
//! including source acceleration and jerk. The original admitted quintic path
//! supplies all three motion derivatives; no trajectory is advanced here.
//! At an internal C2 trajectory knot, jerk and pressure use the left limit;
//! the first knot uses the right limit. Such knots need not be C3 or bandlimited.
//!
//! The medium is homogeneous, quiescent and lossless. No finite-cell quadrature,
//! initial particle-velocity field, absorption or feedback is inferred.
//! Original fixed/moving receiver and ideal-image owners below retain their
//! admitted causal roots and image parity. A moving receiver samples the
//! Eulerian field at its actual position; its observed clock rate changes,
//! but the pressure amplitude receives no receiver multiplier.

use super::TransientQuadrupoleSource;
use crate::AeroacError;
use crate::curle3d::moving::trajectory::PreparedSourceTrajectory;
use crate::curle3d::moving::{
    MovingCompactObserver, PreparedMovingCompactSample, SourceKinematics, UniformSourceMotion,
};
use crate::curle3d::statistics::{invalid, nonfinite};
use fs_math::dd::Dd;

mod paths;

/// Retarded geometry and constant-work pressure coefficients for one source.
///
/// The 18 coefficients multiply physical xx, yy, zz, xy, xz, yz stress entries
/// and their first two derivatives. Preparation/evaluation allocate nothing.
/// Exact source and observer clocks come from the original compact-source
/// propagation owner, with its same subsonic and clearance admission.
#[derive(Clone, Copy, Debug)]
pub struct PreparedMovingQuadrupoleSample {
    propagation: PreparedMovingCompactSample,
    pressure_coefficients: [[f64; 6]; 3],
    direction: [f64; 3],
    sound_speed_m_s: f64,
    maximum_mach: f64,
}

impl MovingCompactObserver {
    /// Prepare complete near/induction/radiation pressure for uniform motion.
    ///
    /// Source values must all be evaluated at the returned emission time.
    /// This does not assign the source a frequency, reconstruct derivatives,
    /// advance its mechanics, or establish finite-history coverage.
    pub fn prepare_uniform_quadrupole(
        &self,
        receiver_time_s: f64,
        motion: &UniformSourceMotion,
        maximum_time_residual_s: f64,
    ) -> Result<PreparedMovingQuadrupoleSample, AeroacError> {
        let propagation = self.prepare_uniform(receiver_time_s, motion, maximum_time_residual_s)?;
        let source = motion.kinematics_at(propagation.source_time_s())?;
        self.prepare_kinematic_quadrupole(propagation, source, [0.0; 3])
    }

    /// Pressure at the arrival of one accelerating source event.
    ///
    /// Position, velocity, acceleration and jerk [m/s³] must describe the same
    /// physical source event in medium-frame axes. The original observer
    /// admits instantaneous speed and clearance; continuous path admission
    /// belongs to the trajectory owner. The tensor history is queried at this
    /// source time, not the returned receiver time.
    pub fn prepare_emission_quadrupole(
        &self,
        source_time_s: f64,
        source: SourceKinematics,
        jerk_m_s3: [f64; 3],
    ) -> Result<PreparedMovingQuadrupoleSample, AeroacError> {
        let propagation = self.prepare_emission(source_time_s, source)?;
        self.prepare_kinematic_quadrupole(propagation, source, jerk_m_s3)
    }

    fn prepare_kinematic_quadrupole(
        &self,
        propagation: PreparedMovingCompactSample,
        source: SourceKinematics,
        jerk_m_s3: [f64; 3],
    ) -> Result<PreparedMovingQuadrupoleSample, AeroacError> {
        if jerk_m_s3.iter().any(|v| !v.is_finite()) {
            return Err(nonfinite("moving quadrupole source jerk"));
        }
        let range = propagation.distance_m();
        let direction: [f64; 3] = core::array::from_fn(|i| {
            (self.observer_position_m()[i] - source.position_m[i]) / range
        });
        let c = self.medium().sound_speed_m_s();
        let mach = source.velocity_m_s.map(|v| v / c);
        let speed = fs_math::det::hypot(fs_math::det::hypot(mach[0], mach[1]), mach[2]);
        let b = (1.0 - speed) * (1.0 + speed);
        let kappa = 1.0 - propagation.radial_mach();
        let uniform = source.acceleration_m_s2 == [0.0; 3] && jerk_m_s3 == [0.0; 3];
        let acceleration = source.acceleration_m_s2.map(|a| (a / c) / c);
        let jerk = jerk_m_s3.map(|j| ((j / c) / c) / c);
        let curvature = b + range * dot(direction, acceleration);
        let jerk_term = range * dot(direction, jerk) - 3.0 * dot(mach, acceleration);
        // D = R kappa, C = 1-|M|² + r.a/c², D_i = C n_i/kappa - M_i.
        // All geometry is evaluated at emission. Differentiating the retarded
        // constraint gives the symmetric tau Hessian below; differentiating
        // D_i adds both acceleration and (r.jerk - 3 v.a)/c³ terms.
        let d_gradient: [f64; 3] =
            core::array::from_fn(|i| curvature * direction[i] / kappa - mach[i]);
        let tau_gradient = direction.map(|n| -(n / c) / kappa);
        let inverse_d = (1.0 / range) / kappa;
        let green = (0.25 / core::f64::consts::PI) * inverse_d;
        let mut coefficients = [[0.0; 6]; 3];
        for (entry, (i, j)) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
            .into_iter()
            .enumerate()
        {
            let (d_hessian, tau_hessian) = if uniform {
                // Preserve the analytic uniform path's original arithmetic.
                let metric = (if i == j { b } else { 0.0 }) + mach[i] * mach[j];
                let d_hessian = (metric - d_gradient[i] * d_gradient[j]) * inverse_d;
                (d_hessian, -(d_hessian / c) / b)
            } else {
                let normal_product = (direction[i] / kappa) * (direction[j] / kappa);
                let metric = (if i == j { 1.0 } else { 0.0 })
                    + (mach[i] * direction[j] + mach[j] * direction[i]) / kappa
                    - curvature * normal_product;
                let tau_hessian = -(metric / c) * inverse_d;
                let d_hessian = curvature * metric * inverse_d
                    + (direction[i] * acceleration[j] + acceleration[i] * direction[j]) / kappa
                    - jerk_term * normal_product;
                (d_hessian, tau_hessian)
            };
            let multiplicity = if i == j { 1.0 } else { 2.0 };
            coefficients[0][entry] = multiplicity
                * (green * inverse_d)
                * (2.0 * d_gradient[i] * d_gradient[j] * inverse_d - d_hessian);
            coefficients[1][entry] = multiplicity
                * green
                * (tau_hessian
                    - inverse_d
                        * (d_gradient[i] * tau_gradient[j] + d_gradient[j] * tau_gradient[i]));
            coefficients[2][entry] = multiplicity * green * tau_gradient[i] * tau_gradient[j];
        }
        if !b.is_finite()
            || b <= 0.0
            || !inverse_d.is_finite()
            || inverse_d <= 0.0
            || !green.is_finite()
            || green <= 0.0
            || coefficients.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "moving quadrupole pressure coefficients are unrepresentable",
            ));
        }
        Ok(PreparedMovingQuadrupoleSample {
            propagation,
            pressure_coefficients: coefficients,
            direction,
            sound_speed_m_s: c,
            maximum_mach: self.limits().maximum_mach,
        })
    }
}

impl PreparedSourceTrajectory {
    /// Observe an admitted accelerating source on the original retarded clock.
    ///
    /// Reuses the trajectory's safeguarded inverse, residual/work limits and
    /// continuous subsonic/clearance admission. Acceleration and analytic jerk
    /// are derivatives of that same quintic, with the documented one-sided
    /// knot convention. Source Q/Qdot/Qddot still come from the caller.
    pub fn prepare_quadrupole_receiver(
        &self,
        receiver_time_s: f64,
        maximum_time_residual_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingQuadrupoleSample, AeroacError> {
        let propagation =
            self.prepare_receiver(receiver_time_s, maximum_time_residual_s, cancelled)?;
        let (source, jerk) = self.kinematics_with_jerk_at(propagation.source_time_s())?;
        self.observer()
            .prepare_kinematic_quadrupole(propagation, source, jerk)
    }
}

impl PreparedMovingQuadrupoleSample {
    /// Original retarded geometry, residual and stationary-receiver clock rate.
    pub const fn propagation(&self) -> &PreparedMovingCompactSample {
        &self.propagation
    }

    /// Physical source-clock query time, before any user history interpolation.
    pub const fn source_time_s(&self) -> f64 {
        self.propagation.source_time_s()
    }

    /// Frequency ratio for a microphone moving in the same fluid rest frame.
    ///
    /// `d tau/d t = (1 - n.v_receiver/c)/(1 - n.v_source/c)`. It describes
    /// timing only and must never multiply the Eulerian pressure below.
    /// Velocity is sampled at reception, source motion at emission. The same
    /// Mach admission as the observer applies to this optional receiver.
    pub fn observed_source_time_rate(
        &self,
        receiver_velocity_m_s: [f64; 3],
    ) -> Result<f64, AeroacError> {
        let mach = receiver_velocity_m_s.map(|v| v / self.sound_speed_m_s);
        if mach.iter().any(|v| !v.is_finite()) {
            return Err(nonfinite("moving quadrupole receiver velocity"));
        }
        let speed = fs_math::det::hypot(fs_math::det::hypot(mach[0], mach[1]), mach[2]);
        if speed > self.maximum_mach {
            return Err(invalid(
                "moving quadrupole receiver exceeds declared Mach limit",
            ));
        }
        let rate = (1.0 - dot(self.direction, mach)) * self.propagation.source_time_rate();
        if !rate.is_finite() || rate <= 0.0 {
            return Err(invalid(
                "moving quadrupole observed clock rate is unrepresentable",
            ));
        }
        Ok(rate)
    }

    /// Evaluate coherent pressure [Pa] with the original physical tensor units.
    ///
    /// `stress_integral_n_m_s` is unused: pressure depends only on Q, Qdot and
    /// Qddot. Particle velocity needs separate causal-history information and
    /// is deliberately not reconstructed as p/Z. Fixed-order DD reduction
    /// retains cancellation between tensor components and radial orders.
    pub fn pressure_pa(&self, source: TransientQuadrupoleSource) -> Result<f64, AeroacError> {
        self.pressure_dd(source).map(Dd::to_f64)
    }

    /// Retain the original compensated pressure through coherent source/image
    /// sums. Round only after the complete receiver event has been assembled.
    /// Source units, finite checks and pressure coefficients match `pressure_pa`.
    pub fn pressure_dd(&self, source: TransientQuadrupoleSource) -> Result<Dd, AeroacError> {
        let values = [
            source.stress_moment_n_m,
            source.stress_rate_n_m_s,
            source.stress_acceleration_n_m_s2,
        ];
        let mut pressure = Dd::ZERO;
        for (coefficients, values) in self.pressure_coefficients.iter().zip(values) {
            for (&coefficient, value) in coefficients.iter().zip(values) {
                if !value.is_finite() {
                    return Err(nonfinite("moving quadrupole stress history"));
                }
                pressure = pressure + Dd::from_f64(coefficient) * Dd::from_f64(value);
            }
        }
        if !pressure.to_f64().is_finite() {
            return Err(invalid("moving quadrupole pressure is unrepresentable"));
        }
        Ok(pressure)
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
