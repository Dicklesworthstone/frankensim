//! Uniformly translating Lighthill stress moments on their original clock.
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
//! Translation is uniform and the medium homogeneous, quiescent and lossless.
//! No acceleration, finite-cell quadrature, initial particle-velocity field,
//! reflections, absorption or source feedback is inferred. A moving receiver
//! may sample the Eulerian field at its actual position; its observed clock
//! rate changes, but the pressure amplitude receives no receiver multiplier.

use super::TransientQuadrupoleSource;
use crate::AeroacError;
use crate::curle3d::moving::{
    MovingCompactObserver, PreparedMovingCompactSample, UniformSourceMotion,
};
use crate::curle3d::statistics::{invalid, nonfinite};
use fs_math::dd::Dd;

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
        let range = propagation.distance_m();
        let direction: [f64; 3] = core::array::from_fn(|i| {
            (self.observer_position_m()[i] - source.position_m[i]) / range
        });
        let c = self.medium().sound_speed_m_s();
        let mach = motion.velocity_m_s().map(|v| v / c);
        let speed = fs_math::det::hypot(fs_math::det::hypot(mach[0], mach[1]), mach[2]);
        let b = (1.0 - speed) * (1.0 + speed);
        let kappa = 1.0 - propagation.radial_mach();
        // D_i = b*n_i/kappa - M_i. Expressing derivatives at the emission
        // event avoids subtracting two large laboratory positions at t.
        let d_gradient: [f64; 3] = core::array::from_fn(|i| b * direction[i] / kappa - mach[i]);
        let tau_gradient = direction.map(|n| -(n / c) / kappa);
        let inverse_d = (1.0 / range) / kappa;
        let green = (0.25 / core::f64::consts::PI) * inverse_d;
        let mut coefficients = [[0.0; 6]; 3];
        for (entry, (i, j)) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
            .into_iter()
            .enumerate()
        {
            let metric = (if i == j { b } else { 0.0 }) + mach[i] * mach[j];
            let d_hessian = (metric - d_gradient[i] * d_gradient[j]) * inverse_d;
            let tau_hessian = -(d_hessian / c) / b;
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
        let pressure = pressure.to_f64();
        if !pressure.is_finite() {
            return Err(invalid("moving quadrupole pressure is unrepresentable"));
        }
        Ok(pressure)
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
