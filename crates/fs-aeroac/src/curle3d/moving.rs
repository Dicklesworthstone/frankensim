//! Retarded pressure of prescribed moving compact 3D sources.
//!
//! In a homogeneous quiescent lossless medium, volume injection Q [m³/s] and
//! force ON the fluid F [N] give
//! `4πp = rho ∂t[Q/(R κ)]ret - div[F/(R κ)]ret`, with
//! `κ = 1 - n.v/c`, `R = |x-y(τ)|`, and `t = τ + R/c`.
//! The Eulerian derivatives retain acceleration, convection and induction;
//! changing only a stationary source's delay would omit these terms. See
//! Farassat, NASA TP-3428, §4.3.1 and NASA/TM-2007-214853:
//! <https://ntrs.nasa.gov/citations/19940029887> and
//! <https://ntrs.nasa.gov/citations/20070010579>.
//!
//! Source mechanics, finite-size compactness, sampled history, antialiasing,
//! scattering, absorption and feedback remain caller-owned. Force pressure
//! has no additional density factor. Transient force particle velocity needs
//! history/initial conditions; this pressure operator does not invent it.

use super::statistics::{invalid, nonfinite};
use crate::AeroacError;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{dd::Dd, det};

/// Uniform-fluid frames and cached fixed laboratory source/microphone paths.
pub mod flow;
/// Finite C2 source trajectories with whole-interval motion/clearance admission.
pub mod trajectory;

/// Physical source and its fixed-frame derivatives at one emission time.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CompactSourceStrength {
    /// Outward volume injection [m³/s], not mass injection.
    pub volume_flow_m3_s: f64,
    /// Time derivative of volume injection [m³/s²].
    pub volume_flow_rate_m3_s2: f64,
    /// Total force ON the fluid [N], not force per unit span.
    pub force_on_fluid_n: [f64; 3],
    /// Fixed-frame time derivative of force ON the fluid [N/s].
    pub force_rate_on_fluid_n_s: [f64; 3],
}

/// A source-history endpoint on its original physical clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TimedCompactSource {
    /// Accepted emission time [s].
    pub source_time_s: f64,
    /// Values and derivatives at exactly this time.
    pub strength: CompactSourceStrength,
}

/// Prescribed source kinematics in the medium's rest frame at emission.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SourceKinematics {
    /// Source position [m].
    pub position_m: [f64; 3],
    /// Source velocity [m/s].
    pub velocity_m_s: [f64; 3],
    /// Source acceleration [m/s²].
    pub acceleration_m_s2: [f64; 3],
}

/// Explicit subsonic margin and geometric clearance, not a compactness test.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingCompactLimits {
    /// Admitted speed/c, inclusive; must satisfy `0 <= maximum_mach < 1`.
    pub maximum_mach: f64,
    /// Minimum source/receiver separation [m], positive and finite.
    pub minimum_distance_m: f64,
}

/// A complete straight-line trajectory, permitting analytic retarded inversion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniformSourceMotion {
    reference_time_s: f64,
    position_m: [f64; 3],
    velocity_m_s: [f64; 3],
}

impl UniformSourceMotion {
    /// Define position and constant velocity on one physical time origin.
    /// Subsonic admission belongs to the observer's original medium.
    pub fn new(
        reference_time_s: f64,
        position_m: [f64; 3],
        velocity_m_s: [f64; 3],
    ) -> Result<Self, AeroacError> {
        if !reference_time_s.is_finite()
            || position_m
                .iter()
                .chain(&velocity_m_s)
                .any(|v| !v.is_finite())
        {
            return Err(nonfinite("uniform source motion"));
        }
        Ok(Self {
            reference_time_s,
            position_m,
            velocity_m_s,
        })
    }

    /// Evaluate the same prescribed trajectory, refusing unrepresentable positions.
    pub fn kinematics_at(&self, time_s: f64) -> Result<SourceKinematics, AeroacError> {
        if !time_s.is_finite() {
            return Err(nonfinite("uniform source time"));
        }
        let dt = Dd::from_f64(time_s) - Dd::from_f64(self.reference_time_s);
        let position_m = core::array::from_fn(|i| {
            (Dd::from_f64(self.position_m[i]) + Dd::from_f64(self.velocity_m_s[i]) * dt).to_f64()
        });
        if position_m.iter().any(|v| !v.is_finite()) {
            return Err(nonfinite("uniform source position"));
        }
        Ok(SourceKinematics {
            position_m,
            velocity_m_s: self.velocity_m_s,
            acceleration_m_s2: [0.0; 3],
        })
    }

    /// Constant velocity [m/s].
    pub const fn velocity_m_s(&self) -> [f64; 3] {
        self.velocity_m_s
    }
    /// Original physical reference time [s].
    pub const fn reference_time_s(&self) -> f64 {
        self.reference_time_s
    }
    /// Position at the reference time [m].
    pub const fn reference_position_m(&self) -> [f64; 3] {
        self.position_m
    }
}

/// One stationary microphone in the original medium's rest frame.
#[derive(Clone, Copy, Debug)]
pub struct MovingCompactObserver {
    medium: FluidAcousticState,
    observer_m: [f64; 3],
    limits: MovingCompactLimits,
}

impl MovingCompactObserver {
    /// Admit the fixed receiver geometry and explicit motion/clearance limits.
    pub fn new(
        medium: FluidAcousticState,
        observer_m: [f64; 3],
        limits: MovingCompactLimits,
    ) -> Result<Self, AeroacError> {
        if observer_m.iter().any(|v| !v.is_finite()) {
            return Err(nonfinite("moving compact receiver position"));
        }
        if !limits.maximum_mach.is_finite()
            || !(0.0..1.0).contains(&limits.maximum_mach)
            || !limits.minimum_distance_m.is_finite()
            || limits.minimum_distance_m <= 0.0
        {
            return Err(invalid("moving compact subsonic margin and clearance"));
        }
        Ok(Self {
            medium,
            observer_m,
            limits,
        })
    }

    /// Original coherent fluid state; no independent rho/c replacement.
    pub const fn medium(&self) -> FluidAcousticState {
        self.medium
    }
    /// Stationary microphone position [m].
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }
    /// Original admission limits.
    pub const fn limits(&self) -> MovingCompactLimits {
        self.limits
    }

    /// Prepare eight pressure coefficients at one supplied emission event.
    ///
    /// Kinematics must all describe this same time. Admission checks the
    /// instantaneous speed only: arbitrary acceleration samples establish no
    /// trajectory-wide bound, uniqueness of retarded roots, or history coverage.
    pub fn prepare_emission(
        &self,
        source_time_s: f64,
        kinematics: SourceKinematics,
    ) -> Result<PreparedMovingCompactSample, AeroacError> {
        if !source_time_s.is_finite()
            || kinematics
                .position_m
                .iter()
                .chain(&kinematics.velocity_m_s)
                .chain(&kinematics.acceleration_m_s2)
                .any(|v| !v.is_finite())
        {
            return Err(nonfinite("moving compact emission"));
        }
        let c = self.medium.sound_speed_m_s();
        let mach = kinematics.velocity_m_s.map(|v| v / c);
        let speed = norm(mach);
        self.admit_speed(speed)?;
        let displacement = core::array::from_fn(|i| self.observer_m[i] - kinematics.position_m[i]);
        let distance_m = norm(displacement);
        if !distance_m.is_finite() || distance_m < self.limits.minimum_distance_m {
            return Err(invalid("moving compact source/receiver separation"));
        }
        let direction = displacement.map(|v| v / distance_m);
        let radial_mach = dot(direction, mach);
        let kappa = 1.0 - radial_mach;
        let source_time_rate = 1.0 / kappa;
        let propagation_delay_s = distance_m / c;
        let g1 = (0.25 / core::f64::consts::PI) / distance_m / kappa / kappa;
        let g2 = g1 / distance_m / kappa;
        let g3 = g1 / kappa;
        if [kappa, source_time_rate, propagation_delay_s, g1, g2, g3]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid("moving compact propagation coefficients"));
        }
        let mach_squared = dot(mach, mach);
        let acceleration_over_c = dot(direction, kinematics.acceleration_m_s2.map(|v| v / c));
        let rho = self.medium.density_kg_m3();
        let volume_rate_scale = rho * g1;
        let force_rate_scale = g1 / c;
        if [volume_rate_scale, force_rate_scale]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid("moving compact source-rate coefficients"));
        }
        // Split radial powers avoid forming R*a or c*c unnecessarily.
        let force_radial = g2 * (1.0 - mach_squared) + (g3 / c) * acceleration_over_c;
        let coefficients = [
            rho * (g3 * acceleration_over_c + (g2 * c) * (radial_mach - mach_squared)),
            volume_rate_scale,
            direction[0] * force_radial - (g1 / distance_m) * mach[0],
            direction[1] * force_radial - (g1 / distance_m) * mach[1],
            direction[2] * force_radial - (g1 / distance_m) * mach[2],
            direction[0] * force_rate_scale,
            direction[1] * force_rate_scale,
            direction[2] * force_rate_scale,
        ];
        let receiver_time_s =
            (Dd::from_f64(source_time_s) + Dd::from_f64(propagation_delay_s)).to_f64();
        if coefficients.iter().any(|v| !v.is_finite())
            || !receiver_time_s.is_finite()
            || receiver_time_s <= source_time_s
        {
            return Err(invalid(
                "moving compact pressure or arrival time is unrepresentable",
            ));
        }
        Ok(PreparedMovingCompactSample {
            coefficients,
            sample: MovingCompactSample {
                source_time_s,
                receiver_time_s,
                propagation_delay_s,
                distance_m,
                radial_mach,
                source_time_rate,
                retarded_time_residual_s: residual(
                    receiver_time_s,
                    source_time_s,
                    propagation_delay_s,
                ),
                monopole_pressure_pa: 0.0,
                dipole_pressure_pa: 0.0,
                pressure_pa: 0.0,
            },
        })
    }

    /// Invert a globally subsonic straight-line trajectory at a receiver time.
    ///
    /// The scaled analytic root avoids squaring physical distances/speeds and
    /// uses the conjugate root for a receding source to prevent cancellation.
    /// `maximum_time_residual_s` is finite and nonnegative; the measured absolute
    /// light-cone residual must fit it. No source-history extrapolation is done
    /// here: the consumer must retain strengths covering the returned source time.
    pub fn prepare_uniform(
        &self,
        receiver_time_s: f64,
        motion: &UniformSourceMotion,
        maximum_time_residual_s: f64,
    ) -> Result<PreparedMovingCompactSample, AeroacError> {
        if !receiver_time_s.is_finite()
            || !maximum_time_residual_s.is_finite()
            || maximum_time_residual_s < 0.0
        {
            return Err(invalid("uniform retarded time and residual budget"));
        }
        let c = self.medium.sound_speed_m_s();
        let mach = motion.velocity_m_s.map(|v| v / c);
        let speed = norm(mach);
        self.admit_speed(speed)?;
        let present = motion.kinematics_at(receiver_time_s)?;
        let displacement = core::array::from_fn(|i| self.observer_m[i] - present.position_m[i]);
        let distance = norm(displacement);
        if !distance.is_finite() || distance <= 0.0 {
            return Err(invalid("uniform retarded displacement"));
        }
        let b = dot(displacement.map(|v| v / distance), mach);
        let a = (1.0 - speed) * (1.0 + speed);
        let discriminant = det::hypot(b, det::sqrt(a));
        let factor = if b < 0.0 {
            1.0 / (discriminant - b)
        } else {
            (discriminant + b) / a
        };
        let delay = (distance / c) * factor;
        let source_time_s = (Dd::from_f64(receiver_time_s) - Dd::from_f64(delay)).to_f64();
        if !delay.is_finite()
            || delay <= 0.0
            || !source_time_s.is_finite()
            || source_time_s >= receiver_time_s
        {
            return Err(invalid("uniform retarded travel time is unrepresentable"));
        }
        let mut prepared =
            self.prepare_emission(source_time_s, motion.kinematics_at(source_time_s)?)?;
        let error = residual(
            receiver_time_s,
            source_time_s,
            prepared.sample.propagation_delay_s,
        );
        if !error.is_finite() || error > maximum_time_residual_s {
            return Err(AeroacError::NotConverged {
                what: "uniform retarded-time residual",
                residual: error,
            });
        }
        prepared.sample.receiver_time_s = receiver_time_s;
        prepared.sample.retarded_time_residual_s = error;
        Ok(prepared)
    }

    fn admit_speed(&self, speed: f64) -> Result<(), AeroacError> {
        if !speed.is_finite() || speed > self.limits.maximum_mach {
            Err(invalid("moving compact source exceeds declared Mach limit"))
        } else {
            Ok(())
        }
    }
}

/// Prepared moving geometry, with no allocation or history owned by this value.
#[derive(Clone, Copy, Debug)]
pub struct PreparedMovingCompactSample {
    coefficients: [f64; 8],
    sample: MovingCompactSample,
}

impl PreparedMovingCompactSample {
    /// Emission time at which the caller must evaluate all source quantities [s].
    pub const fn source_time_s(&self) -> f64 {
        self.sample.source_time_s
    }
    /// Arrival time, or the exact requested receiver time for uniform inversion [s].
    pub const fn receiver_time_s(&self) -> f64 {
        self.sample.receiver_time_s
    }
    /// Physical geometric travel time R/c [s].
    pub const fn propagation_delay_s(&self) -> f64 {
        self.sample.propagation_delay_s
    }
    /// Emission-to-receiver separation [m].
    pub const fn distance_m(&self) -> f64 {
        self.sample.distance_m
    }
    /// Emission velocity projected toward the receiver, divided by sound speed.
    pub const fn radial_mach(&self) -> f64 {
        self.sample.radial_mach
    }
    /// Local emission-time rate dτ/dt = 1/(1-Mr) at the stationary receiver.
    pub const fn source_time_rate(&self) -> f64 {
        self.sample.source_time_rate
    }
    /// Absolute residual |t-τ-R/c| [s], evaluated with double-double subtraction.
    pub const fn retarded_time_residual_s(&self) -> f64 {
        self.sample.retarded_time_residual_s
    }

    /// Apply coherent strengths at exactly the prepared emission time.
    /// Fixed-order double-double products/sums preserve cancellation residuals;
    /// all eight components are bounded constant work, without allocation.
    pub fn evaluate(
        &self,
        strength: CompactSourceStrength,
    ) -> Result<MovingCompactSample, AeroacError> {
        let (sample, _) = self.evaluate_with_pressure_dd(strength)?;
        Ok(sample)
    }

    // Keep weak source components available to the coherent image-path sum.
    // The public sample and all per-path overflow checks remain unchanged.
    fn evaluate_with_pressure_dd(
        &self,
        strength: CompactSourceStrength,
    ) -> Result<(MovingCompactSample, Dd), AeroacError> {
        let values = [
            strength.volume_flow_m3_s,
            strength.volume_flow_rate_m3_s2,
            strength.force_on_fluid_n[0],
            strength.force_on_fluid_n[1],
            strength.force_on_fluid_n[2],
            strength.force_rate_on_fluid_n_s[0],
            strength.force_rate_on_fluid_n_s[1],
            strength.force_rate_on_fluid_n_s[2],
        ];
        if values.iter().any(|v| !v.is_finite()) {
            return Err(nonfinite("moving compact source strength"));
        }
        let mut total = Dd::ZERO;
        let mut monopole = Dd::ZERO;
        let mut dipole = Dd::ZERO;
        for (i, (coefficient, value)) in self.coefficients.iter().zip(values).enumerate() {
            let term = Dd::from_f64(*coefficient) * Dd::from_f64(value);
            total = total + term;
            if i < 2 {
                monopole = monopole + term;
            } else {
                dipole = dipole + term;
            }
        }
        let mut sample = self.sample;
        sample.monopole_pressure_pa = monopole.to_f64();
        sample.dipole_pressure_pa = dipole.to_f64();
        sample.pressure_pa = total.to_f64();
        if [
            sample.monopole_pressure_pa,
            sample.dipole_pressure_pa,
            sample.pressure_pa,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(nonfinite("moving compact pressure"));
        }
        Ok((sample, total))
    }
}

/// Physical pressure with the actual emission/receiver event and its residual.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingCompactSample {
    /// Source evaluation time [s].
    pub source_time_s: f64,
    /// Receiver evaluation time [s].
    pub receiver_time_s: f64,
    /// Geometric propagation delay R/c [s].
    pub propagation_delay_s: f64,
    /// Geometric separation at emission [m].
    pub distance_m: f64,
    /// Radial emission Mach number n.v/c.
    pub radial_mach: f64,
    /// Local dτ/dt = 1/(1-Mr), including Doppler time compression.
    pub source_time_rate: f64,
    /// Absolute light-cone residual |t-τ-R/c| [s].
    pub retarded_time_residual_s: f64,
    /// Signed volume-source pressure [Pa].
    pub monopole_pressure_pa: f64,
    /// Signed force-source pressure [Pa], including induction.
    pub dipole_pressure_pa: f64,
    /// Coherent total pressure [Pa], before audio gain or normalization.
    pub pressure_pa: f64,
}

fn norm(vector: [f64; 3]) -> f64 {
    det::hypot(det::hypot(vector[0], vector[1]), vector[2])
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    let mut sum = Dd::ZERO;
    for i in 0..3 {
        sum = sum + Dd::from_f64(a[i]) * Dd::from_f64(b[i]);
    }
    sum.to_f64()
}

fn residual(receiver: f64, source: f64, delay: f64) -> f64 {
    (Dd::from_f64(receiver) - Dd::from_f64(source) - Dd::from_f64(delay))
        .abs()
        .to_f64()
}
