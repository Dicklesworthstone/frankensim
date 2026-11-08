//! Complete moving-stress fields with an explicit past-worldline integral.
//!
//! Let H_ij(x,t) = integral_{t0}^{tau(x,t)} Q_ij(s)/|x-y(s)| ds.
//! Away from the complete source worldline, linear momentum gives
//! u_a = -partial_a partial_i partial_j H_ij / (4 pi rho).
//! The original pressure remains the Hessian of the retarded Green potential.
//! The source analogy and retarded Jacobian are Lighthill (1952),
//! <https://doi.org/10.1098/rspa.1952.0060>, and Farassat (2007),
//! <https://ntrs.nasa.gov/citations/20070010579>.
//!
//! Differentiating the moving upper limit three times gives local terms in
//! Q, Qdot and Qddot, and the nonlocal integral of
//! [15 n (n.Q.n) - 6 Q n - 3 tr(Q) n]/R^4.
//! Every history point uses the SAME current receiver position. The six
//! scalar integrals of Q alone cannot supply this moving-source history.
//!
//! Local coefficients below are direct Leibniz/implicit derivatives, not
//! numerical differences. The dimensionless clock derivatives are c*tau_i,
//! c*R*tau_ij and c*R^2*tau_aij. Their analytic construction includes source
//! acceleration and jerk and requires no higher stress or motion derivative.
//! The original pressure-only preparation and arithmetic remain unchanged.

use super::{PreparedMovingQuadrupoleSample, dot};
use crate::AeroacError;
use crate::curle3d::moving::trajectory::PreparedSourceTrajectory;
use crate::curle3d::moving::{
    MovingCompactObserver, PreparedMovingCompactSample, SourceKinematics, UniformSourceMotion,
};
use crate::curle3d::quadrupole::transient::{RetardedQuadrupoleField, TransientQuadrupoleSource};
use crate::curle3d::statistics::{invalid, nonfinite};
use fs_math::{dd::Dd, det};

/// Physical authority retained by a full moving-stress field.
pub const MOVING_QUADRUPOLE_FIELD_SCOPE: &str = "Prescribed moving symmetric 3D Lighthill stress moments [N m] on their original clock; complete pressure and particle velocity in a homogeneous lossless fluid rest frame, with explicit past-worldline stress integral and original ideal images when supplied. No inferred initial field, finite-cell quadrature, unresolved turbulence, absorption, general scattering, source feedback or measured-SPL qualification.";

const ENTRIES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

/// Constant-work local field coefficients for one original retarded event.
///
/// The separate history integral must be evaluated at this exact receiver
/// position and integrated on the original source clock through this emission
/// event. Zero history at a declared origin means zero prior stress forcing;
/// it does not supply an unknown homogeneous initial acoustic field.
///
/// Ideal-image preparations consume an UNTRANSFORMED history integral at
/// their virtual receiver, using the ORIGINAL source tensor/worldline. Their
/// cached velocity map applies parity and frame transformation exactly once.
/// No new storage or work is added to the existing pressure-only sample.
#[derive(Clone, Copy, Debug)]
pub struct PreparedMovingQuadrupoleFieldSample {
    pressure: PreparedMovingQuadrupoleSample,
    // Cartesian velocity, stress derivative order, symmetric tensor entry.
    velocity_coefficients: [[[f64; 6]; 3]; 3],
    history_velocity: [[f64; 3]; 3],
}

impl MovingCompactObserver {
    // A coherent Q/F/stress composition already owns this admitted causal root
    // and its original motion jet. Keep both tensor derivatives in this owner
    // without another retarded search or reconstructing motion at a new time.
    pub(in crate::curle3d) fn prepare_quadrupole_field_sample(
        &self,
        propagation: PreparedMovingCompactSample,
        source: SourceKinematics,
        jerk_m_s3: [f64; 3],
    ) -> Result<PreparedMovingQuadrupoleFieldSample, AeroacError> {
        let pressure = self.prepare_kinematic_quadrupole(propagation, source, jerk_m_s3)?;
        self.prepare_kinematic_quadrupole_field(pressure, source, jerk_m_s3)
    }

    /// Full-field coefficients for uniform motion on the original causal root.
    pub fn prepare_uniform_quadrupole_field(
        &self,
        receiver_time_s: f64,
        motion: &UniformSourceMotion,
        maximum_time_residual_s: f64,
    ) -> Result<PreparedMovingQuadrupoleFieldSample, AeroacError> {
        let propagation = self.prepare_uniform(receiver_time_s, motion, maximum_time_residual_s)?;
        let source = motion.kinematics_at(propagation.source_time_s())?;
        let pressure = self.prepare_kinematic_quadrupole(propagation, source, [0.0; 3])?;
        self.prepare_kinematic_quadrupole_field(pressure, source, [0.0; 3])
    }

    /// Full-field coefficients at the arrival of one accelerating source event.
    ///
    /// All motion derivatives belong to this same event in medium-frame axes.
    /// Past-worldline clearance and causal stress history are separate inputs.
    pub fn prepare_emission_quadrupole_field(
        &self,
        source_time_s: f64,
        source: SourceKinematics,
        jerk_m_s3: [f64; 3],
    ) -> Result<PreparedMovingQuadrupoleFieldSample, AeroacError> {
        let pressure = self.prepare_emission_quadrupole(source_time_s, source, jerk_m_s3)?;
        self.prepare_kinematic_quadrupole_field(pressure, source, jerk_m_s3)
    }

    // Native motion owners pass the already solved/admitted root and exact
    // original kinematics. There is no second inverse or pressure evaluation.
    pub(super) fn prepare_kinematic_quadrupole_field(
        &self,
        pressure: PreparedMovingQuadrupoleSample,
        source: SourceKinematics,
        jerk_m_s3: [f64; 3],
    ) -> Result<PreparedMovingQuadrupoleFieldSample, AeroacError> {
        let radius = pressure.propagation().distance_m();
        let n = pressure.direction;
        let c = self.medium().sound_speed_m_s();
        let rho = self.medium().density_kg_m3();
        let mach = source.velocity_m_s.map(|v| v / c);
        let acceleration = source.acceleration_m_s2.map(|a| ((a / c) / c) * radius);
        let jerk = jerk_m_s3.map(|j| (((j / c) / c) / c) * radius * radius);
        let radial_mach = pressure.propagation().radial_mach();
        let kappa = 1.0 - radial_mach;
        let history_scale = (0.25 / core::f64::consts::PI) / rho;
        // Split inverse powers; never form rho*c or a large R^3.
        let scales = [
            history_scale / radius / radius / radius / c,
            history_scale / radius / radius / c / c,
            history_scale / radius / c / c / c,
        ];
        if acceleration.iter().chain(&jerk).any(|v| !v.is_finite())
            || !history_scale.is_finite()
            || history_scale <= 0.0
            || scales.iter().any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "moving quadrupole velocity geometry is unrepresentable",
            ));
        }
        let clock = RetardedClock::new(n, mach, acceleration, jerk, radial_mach, kappa);
        let green_time = radial_mach;
        let green_space_time: [f64; 3] =
            core::array::from_fn(|i| mach[i] - 3.0 * n[i] * radial_mach);
        let green_time_time =
            3.0 * radial_mach * radial_mach - dot(mach, mach) + dot(n, acceleration);
        let mut velocity_coefficients = [[[0.0; 6]; 3]; 3];
        for (a, orders) in velocity_coefficients.iter_mut().enumerate() {
            for (entry, (i, j)) in ENTRIES.into_iter().enumerate() {
                let t = clock.first;
                let tt = clock.second;
                let product = t[a] * t[i] * t[j];
                let mixed = tt[a][i] * t[j] + tt[a][j] * t[i] + tt[i][j] * t[a];
                let hessian = |i: usize, j: usize| {
                    (3.0 * n[i]).mul_add(n[j], if i == j { -1.0 } else { 0.0 })
                };
                // Third spatial Leibniz derivative of the time integral:
                // integral h_aij plus the differentiated retarded endpoint.
                let local_moment = hessian(i, j) * t[a]
                    + hessian(a, i) * t[j]
                    + hessian(a, j) * t[i]
                    + green_space_time[i] * t[j] * t[a]
                    + green_space_time[j] * t[i] * t[a]
                    + green_space_time[a] * t[i] * t[j]
                    + green_time_time * product
                    - n[i] * tt[j][a]
                    - n[j] * tt[i][a]
                    - n[a] * tt[i][j]
                    + green_time * mixed
                    + clock.third[a][i][j];
                let local_rate = -n[i] * t[j] * t[a] - n[j] * t[i] * t[a] - n[a] * t[i] * t[j]
                    + 2.0 * green_time * product
                    + mixed;
                let multiplicity = if i == j { 1.0 } else { 2.0 };
                for (order, value) in [local_moment, local_rate, product].into_iter().enumerate() {
                    orders[order][entry] = -multiplicity * value * scales[order];
                }
            }
        }
        if velocity_coefficients
            .iter()
            .flatten()
            .flatten()
            .any(|v| !v.is_finite())
        {
            return Err(nonfinite("moving quadrupole local velocity coefficients"));
        }
        Ok(PreparedMovingQuadrupoleFieldSample {
            pressure,
            velocity_coefficients,
            history_velocity: core::array::from_fn(|i| {
                core::array::from_fn(|j| if i == j { history_scale } else { 0.0 })
            }),
        })
    }

    /// Spatially weighted past-stress integrand at this CURRENT receiver.
    ///
    /// Returns [15 n (n.Q.n) - 6 Q n - 3 tr(Q) n]/R^4 [N/m^3]. Integrate
    /// from the declared source-history origin through the retarded emission.
    /// The six entries are xx, yy, zz, xy, xz, yz [N m], with physical shear.
    /// No Doppler factor enters this integrand: the changing upper limit is
    /// already differentiated in the cached local field coefficients.
    ///
    /// This point query checks finite values and separation. The history owner
    /// must separately establish continuous past-worldline clearance between
    /// quadrature nodes and retain the same source clock/receiver/image path.
    pub fn quadrupole_history_integrand(
        &self,
        source_position_m: [f64; 3],
        stress_moment_n_m: [f64; 6],
    ) -> Result<[f64; 3], AeroacError> {
        if source_position_m
            .iter()
            .chain(&stress_moment_n_m)
            .any(|v| !v.is_finite())
        {
            return Err(nonfinite("moving quadrupole history geometry or strength"));
        }
        let displacement: [f64; 3] =
            core::array::from_fn(|i| self.observer_position_m()[i] - source_position_m[i]);
        let radius = det::hypot(
            det::hypot(displacement[0], displacement[1]),
            displacement[2],
        );
        if !radius.is_finite() || radius < self.limits().minimum_distance_m {
            return Err(invalid(
                "moving quadrupole history source/receiver separation",
            ));
        }
        let inverse_fourth = 1.0 / radius / radius / radius / radius;
        if !inverse_fourth.is_finite() || inverse_fourth <= 0.0 {
            return Err(invalid(
                "moving quadrupole history geometry is unrepresentable",
            ));
        }
        let n = displacement.map(|v| Dd::from_f64(v / radius));
        let q = stress_moment_n_m.map(Dd::from_f64);
        // S = 3 Q - tr(Q) I. Differences remove an isotropic tensor exactly,
        // without relying on a rounded direction to have exactly unit norm.
        // The kernel is then [5 n (n.S.n) - 2 S n]/R^4.
        let xx = (q[0] - q[1]) + (q[0] - q[2]);
        let yy = (q[1] - q[0]) + (q[1] - q[2]);
        let s = [
            xx,
            yy,
            -(xx + yy),
            q[3] * Dd::from_f64(3.0),
            q[4] * Dd::from_f64(3.0),
            q[5] * Dd::from_f64(3.0),
        ];
        let tensor = [[s[0], s[3], s[4]], [s[3], s[1], s[5]], [s[4], s[5], s[2]]];
        let sn: [Dd; 3] =
            core::array::from_fn(|i| (0..3).fold(Dd::ZERO, |sum, j| sum + tensor[i][j] * n[j]));
        let axial = (0..3).fold(Dd::ZERO, |sum, i| sum + sn[i] * n[i]);
        let result = core::array::from_fn(|i| {
            ((Dd::from_f64(5.0) * n[i] * axial - Dd::from_f64(2.0) * sn[i])
                * Dd::from_f64(inverse_fourth))
            .to_f64()
        });
        if result.iter().any(|v| !v.is_finite()) {
            return Err(nonfinite("moving quadrupole history integrand"));
        }
        Ok(result)
    }
}

impl PreparedSourceTrajectory {
    /// Full field using the original trajectory's root and analytic motion jet.
    pub fn prepare_quadrupole_field_receiver(
        &self,
        receiver_time_s: f64,
        maximum_time_residual_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingQuadrupoleFieldSample, AeroacError> {
        let propagation =
            self.prepare_receiver(receiver_time_s, maximum_time_residual_s, cancelled)?;
        let (source, jerk) = self.kinematics_with_jerk_at(propagation.source_time_s())?;
        let observer = self.observer();
        let pressure = observer.prepare_kinematic_quadrupole(propagation, source, jerk)?;
        observer.prepare_kinematic_quadrupole_field(pressure, source, jerk)
    }
}

impl PreparedMovingQuadrupoleFieldSample {
    /// The unchanged original pressure operator and its clock-rate query.
    pub const fn pressure(&self) -> &PreparedMovingQuadrupoleSample {
        &self.pressure
    }

    /// Original admitted retarded root, residual and geometry.
    pub const fn propagation(&self) -> &PreparedMovingCompactSample {
        self.pressure.propagation()
    }

    /// Original physical source-clock query time.
    pub const fn source_time_s(&self) -> f64 {
        self.pressure.source_time_s()
    }

    /// Full field, retaining the caller's explicit initial stress contribution.
    ///
    /// `source.stress_integral_n_m_s` is unused: the supplied spatially weighted
    /// worldline integral is required. Sum coherent p/u from different sources
    /// and paths before forming intensity; `evaluate_dd` supports that use.
    pub fn evaluate(
        &self,
        source: TransientQuadrupoleSource,
        history_integral_n_s_m3: [f64; 3],
    ) -> Result<RetardedQuadrupoleField, AeroacError> {
        let (pressure, velocity) =
            self.evaluate_dd(source, history_integral_n_s_m3.map(Dd::from_f64))?;
        let particle_velocity_m_s = velocity.map(Dd::to_f64);
        let instantaneous_intensity_w_m2 = velocity.map(|u| (pressure * u).to_f64());
        if instantaneous_intensity_w_m2.iter().any(|v| !v.is_finite()) {
            return Err(nonfinite("moving quadrupole instantaneous intensity"));
        }
        Ok(RetardedQuadrupoleField {
            pressure_pa: pressure.to_f64(),
            particle_velocity_m_s,
            instantaneous_intensity_w_m2,
            scope: MOVING_QUADRUPOLE_FIELD_SCOPE,
        })
    }

    /// Unrounded first-order fields for coherent source/image accumulation.
    ///
    /// History is an integral at this exact current virtual receiver, without
    /// image parity or a vector transform. The prepared field applies that map
    /// once. Missing past source history is never supplied by this operation.
    pub fn evaluate_dd(
        &self,
        source: TransientQuadrupoleSource,
        history_integral_n_s_m3: [Dd; 3],
    ) -> Result<(Dd, [Dd; 3]), AeroacError> {
        if history_integral_n_s_m3
            .iter()
            .any(|v| !v.to_f64().is_finite())
        {
            return Err(nonfinite("moving quadrupole history integral"));
        }
        let pressure = self.pressure.pressure_dd(source)?;
        let strengths = [
            source.stress_moment_n_m,
            source.stress_rate_n_m_s,
            source.stress_acceleration_n_m_s2,
        ];
        let velocity: [Dd; 3] = core::array::from_fn(|axis| {
            let mut value = Dd::ZERO;
            for (coefficients, values) in self.velocity_coefficients[axis].iter().zip(strengths) {
                for (&coefficient, source) in coefficients.iter().zip(values) {
                    value = value + Dd::from_f64(coefficient) * Dd::from_f64(source);
                }
            }
            for (coefficient, history) in self.history_velocity[axis]
                .into_iter()
                .zip(history_integral_n_s_m3)
            {
                value = value + Dd::from_f64(coefficient) * history;
            }
            value
        });
        if velocity.iter().any(|v| !v.to_f64().is_finite()) {
            return Err(nonfinite("moving quadrupole particle velocity"));
        }
        Ok((pressure, velocity))
    }

    pub(super) fn with_image(
        mut self,
        parity: f64,
        velocity_map: [[f64; 3]; 3],
    ) -> Result<Self, AeroacError> {
        self.pressure = self.pressure.with_image(parity, velocity_map);
        let original = self.velocity_coefficients;
        self.velocity_coefficients = core::array::from_fn(|i| {
            core::array::from_fn(|order| {
                core::array::from_fn(|entry| {
                    (0..3)
                        .fold(Dd::ZERO, |sum, j| {
                            sum + Dd::from_f64(velocity_map[i][j])
                                * Dd::from_f64(original[j][order][entry])
                        })
                        .to_f64()
                })
            })
        });
        let original_history = self.history_velocity;
        self.history_velocity = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (0..3)
                    .fold(Dd::ZERO, |sum, k| {
                        sum + Dd::from_f64(velocity_map[i][k])
                            * Dd::from_f64(original_history[k][j])
                    })
                    .to_f64()
            })
        });
        if self
            .velocity_coefficients
            .iter()
            .flatten()
            .flatten()
            .chain(self.history_velocity.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(nonfinite("moving quadrupole image velocity coefficients"));
        }
        Ok(self)
    }
}

// Implicitly differentiate tau + R(x,tau)/c = t. Scaled source derivatives
// A=R*a/c^2 and J=R^2*jerk/c^3 keep the clock calculation dimensionless.
struct RetardedClock {
    first: [f64; 3],
    second: [[f64; 3]; 3],
    third: [[[f64; 3]; 3]; 3],
}

impl RetardedClock {
    fn new(n: [f64; 3], m: [f64; 3], a: [f64; 3], jerk: [f64; 3], mu: f64, kappa: f64) -> Self {
        let m2 = dot(m, m);
        let na = dot(n, a);
        // Scaled fixed-source-time R derivatives: R_it, R_tt, R_ijt,
        // R_itt and R_ttt. These are NOT derivatives along the retarded root.
        let space_time: [f64; 3] = core::array::from_fn(|i| -m[i] + n[i] * mu);
        let time_time = m2 - mu * mu - na;
        let space_space_time: [[f64; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                m[i] * n[j]
                    + n[i] * m[j]
                    + ((if i == j { 1.0 } else { 0.0 }) - 3.0 * n[i] * n[j]) * mu
            })
        });
        let space_time_time: [f64; 3] = core::array::from_fn(|i| {
            -2.0 * mu * m[i] + (3.0 * mu * mu - m2) * n[i] - a[i] + n[i] * na
        });
        let time_time_time = 3.0 * (dot(m, a) - mu * na) + 3.0 * mu * (m2 - mu * mu) - dot(n, jerk);
        let first = n.map(|v| -v / kappa);
        let second: [[f64; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                -((if i == j { 1.0 } else { 0.0 }) - n[i] * n[j]
                    + space_time[i] * first[j]
                    + space_time[j] * first[i]
                    + time_time * first[i] * first[j])
                    / kappa
            })
        });
        let third = core::array::from_fn(|a| {
            core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    let space = 3.0 * n[a] * n[i] * n[j]
                        - (if a == i { n[j] } else { 0.0 })
                        - (if a == j { n[i] } else { 0.0 })
                        - (if i == j { n[a] } else { 0.0 });
                    -(space
                        + space_space_time[a][i] * first[j]
                        + space_space_time[a][j] * first[i]
                        + space_space_time[i][j] * first[a]
                        + space_time_time[a] * first[i] * first[j]
                        + space_time_time[i] * first[a] * first[j]
                        + space_time_time[j] * first[a] * first[i]
                        + time_time_time * first[a] * first[i] * first[j]
                        + space_time[a] * second[i][j]
                        + space_time[i] * second[a][j]
                        + space_time[j] * second[a][i]
                        + time_time
                            * (second[a][i] * first[j]
                                + second[a][j] * first[i]
                                + second[i][j] * first[a]))
                        / kappa
                })
            })
        });
        Self {
            first,
            second,
            third,
        }
    }
}
