//! Compact-source pressure in a uniform, subsonic, moving fluid.
//!
//! Source and microphone positions are fixed in the laboratory. The original
//! fluid's sound speed is relative to that fluid, not the laboratory. In the
//! fluid rest frame both positions move with velocity `-U`; the existing
//! moving-source pressure owner therefore supplies the complete coefficients,
//! including induction and the convective volume-source term.
//!
//! For `r = x-y`, `M = U/c`, `beta² = 1-|M|²`, define
//! `R* = sqrt(beta² |r|² + (M.r)²)` and `T = (R* - M.r)/(c beta²)`.
//! The convected Green function is `G = exp(i omega T)/(4 pi R*)` under
//! `exp(-i omega t)`. The underlying source equation is
//! `(D²/c² - Laplacian) p = rho D[Q delta_y] - div[F delta_y]`,
//! where `D = d/dt + U.grad`. Thus the harmonic pressure is
//! `rho (-i omega + U.grad) G Q - F.grad(G)`, not a stationary Green
//! function with an adjusted scalar sound speed. Q is volume injection and F
//! is force on the fluid in this linearized momentum equation; any additional
//! momentum carried by an injected jet must be supplied by the source model.
//!
//! The Galilean construction uses the retarded potential differentiated in
//! Farassat, NASA/TM-2007-214853, <https://ntrs.nasa.gov/citations/20070010579>.
//! Both laboratory endpoints are stationary, so the observed frequency equals
//! the source frequency even though propagation is direction dependent.
//! This is a homogeneous lossless free-field pressure operator. Shear/refraction,
//! boundaries, molecular relaxation, particle velocity, acoustic energy flux
//! and source feedback are not supplied by this pressure owner. Moving
//! laboratory paths use [`UniformFlowFrame`] with the original accepted
//! trajectory and source/receiver owners.

use super::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, PreparedMovingCompactSample,
    SourceKinematics, dot, norm, residual,
};
use crate::AeroacError;
use crate::curle3d::CompactSourcePhasor;
use crate::curle3d::statistics::{finite, invalid, nonfinite};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd, det};

/// A shared fluid-rest coordinate frame for accepted laboratory trajectories.
///
/// The Galilean map is `x' = x - origin - U*(t-epoch)`, `v' = v-U`, `a' = a`.
/// Time, volume injection, force vectors and their time derivatives retain
/// their physical units and values. The explicit reference event keeps large
/// absolute laboratory positions and clock epochs out of the interpolant.
/// Both paths of a causal observation must retain the SAME frame.
///
/// This is a constant-velocity transformation, not a spatially varying-flow
/// model. The original trajectory owner admits source and microphone speeds
/// relative to the fluid; changing frames does not relax those Mach limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniformFlowFrame {
    velocity_m_s: [f64; 3],
    reference_time_s: f64,
    reference_position_m: [f64; 3],
}

impl UniformFlowFrame {
    /// Bind a finite fluid velocity and one explicit laboratory reference event.
    pub fn new(
        velocity_m_s: [f64; 3],
        reference_time_s: f64,
        reference_position_m: [f64; 3],
    ) -> Result<Self, AeroacError> {
        if !reference_time_s.is_finite()
            || velocity_m_s
                .iter()
                .chain(&reference_position_m)
                .any(|v| !v.is_finite())
        {
            return Err(nonfinite("uniform-flow reference frame"));
        }
        Ok(Self {
            velocity_m_s,
            reference_time_s,
            reference_position_m,
        })
    }

    /// Uniform fluid velocity relative to the laboratory [m/s].
    pub const fn velocity_m_s(self) -> [f64; 3] {
        self.velocity_m_s
    }
    /// Original laboratory clock epoch [s].
    pub const fn reference_time_s(self) -> f64 {
        self.reference_time_s
    }
    /// Laboratory position assigned zero fluid-frame position at the epoch [m].
    pub const fn reference_position_m(self) -> [f64; 3] {
        self.reference_position_m
    }

    /// Transform one laboratory event before the original quintic compilation.
    pub fn to_fluid_frame(
        self,
        time_s: f64,
        laboratory: SourceKinematics,
    ) -> Result<SourceKinematics, AeroacError> {
        self.transform(time_s, laboratory, false)
    }

    /// Read one interpolated event back in the original laboratory frame.
    /// This observation changes neither accepted knots nor acoustic coefficients.
    pub fn to_laboratory_frame(
        self,
        time_s: f64,
        fluid: SourceKinematics,
    ) -> Result<SourceKinematics, AeroacError> {
        self.transform(time_s, fluid, true)
    }

    fn transform(
        self,
        time_s: f64,
        original: SourceKinematics,
        to_laboratory: bool,
    ) -> Result<SourceKinematics, AeroacError> {
        if !time_s.is_finite()
            || original
                .position_m
                .iter()
                .chain(&original.velocity_m_s)
                .chain(&original.acceleration_m_s2)
                .any(|v| !v.is_finite())
        {
            return Err(nonfinite("uniform-flow kinematic event"));
        }
        let dt = Dd::from_f64(time_s) - Dd::from_f64(self.reference_time_s);
        let mut transformed = original;
        for axis in 0..3 {
            let velocity = Dd::from_f64(self.velocity_m_s[axis]);
            let drift = velocity * dt;
            let origin = Dd::from_f64(self.reference_position_m[axis]);
            let position = Dd::from_f64(original.position_m[axis]);
            let speed = Dd::from_f64(original.velocity_m_s[axis]);
            let (position, speed) = if to_laboratory {
                (position + origin + drift, speed + velocity)
            } else {
                (position - origin - drift, speed - velocity)
            };
            transformed.position_m[axis] = position.to_f64();
            transformed.velocity_m_s[axis] = speed.to_f64();
        }
        if transformed
            .position_m
            .iter()
            .chain(&transformed.velocity_m_s)
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "uniform-flow kinematic transform is unrepresentable",
            ));
        }
        Ok(transformed)
    }
}

/// One fixed laboratory microphone in a homogeneous subsonic fluid flow.
/// Sources bound to this receiver share its original fluid, velocity and frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniformFlowReceiver {
    medium: FluidAcousticState,
    velocity_m_s: [f64; 3],
    observer_m: [f64; 3],
    limits: MovingCompactLimits,
}

impl UniformFlowReceiver {
    /// Admit the physical receiver, fluid speed and laboratory clearance policy.
    pub fn new(
        medium: FluidAcousticState,
        velocity_m_s: [f64; 3],
        observer_m: [f64; 3],
        limits: MovingCompactLimits,
    ) -> Result<Self, AeroacError> {
        let admission = MovingCompactObserver::new(medium, observer_m, limits)?;
        if velocity_m_s.iter().any(|v| !v.is_finite()) {
            return Err(nonfinite("uniform-flow fluid velocity"));
        }
        admission.admit_speed(norm(velocity_m_s.map(|v| v / medium.sound_speed_m_s())))?;
        Ok(Self {
            medium,
            velocity_m_s,
            observer_m,
            limits,
        })
    }

    /// Prepare a source path without changing the common receiver or fluid.
    pub fn bind_source(
        self,
        source_m: [f64; 3],
    ) -> Result<UniformFlowCompactObserver, AeroacError> {
        UniformFlowCompactObserver::from_receiver(self, source_m)
    }
    /// Original coherent thermodynamic acoustic state.
    pub const fn medium(self) -> FluidAcousticState {
        self.medium
    }
    /// Fluid velocity relative to the laboratory [m/s].
    pub const fn velocity_m_s(self) -> [f64; 3] {
        self.velocity_m_s
    }
    /// Fixed microphone position in the laboratory [m].
    pub const fn observer_position_m(self) -> [f64; 3] {
        self.observer_m
    }
    /// Original fluid-Mach and laboratory-clearance limits.
    pub const fn limits(self) -> MovingCompactLimits {
        self.limits
    }
}

/// Cached pressure propagation between fixed laboratory positions in one flow.
///
/// Preparation and every observation are constant-work and allocation-free.
/// All eight real pressure coefficients come from the existing moving-source
/// owner; neither a second acoustic source law nor an effective sound speed is
/// introduced. The physical density and sound speed remain one original state.
#[derive(Clone, Copy, Debug)]
pub struct UniformFlowCompactObserver {
    receiver: UniformFlowReceiver,
    source_m: [f64; 3],
    pressure: PreparedMovingCompactSample,
}

impl UniformFlowCompactObserver {
    /// Admit a homogeneous subsonic flow and fixed laboratory geometry.
    ///
    /// `limits.maximum_mach` bounds the fluid speed/c; `minimum_distance_m`
    /// bounds the actual laboratory source/microphone separation. These are
    /// geometric limits, not finite-source compactness or linearity evidence.
    pub fn new(
        medium: FluidAcousticState,
        velocity_m_s: [f64; 3],
        observer_m: [f64; 3],
        source_m: [f64; 3],
        limits: MovingCompactLimits,
    ) -> Result<Self, AeroacError> {
        UniformFlowReceiver::new(medium, velocity_m_s, observer_m, limits)?.bind_source(source_m)
    }

    fn from_receiver(
        receiver: UniformFlowReceiver,
        source_m: [f64; 3],
    ) -> Result<Self, AeroacError> {
        if source_m.iter().any(|v| !v.is_finite()) {
            return Err(nonfinite("uniform-flow source"));
        }
        let UniformFlowReceiver {
            medium,
            velocity_m_s,
            observer_m,
            limits,
        } = receiver;
        let c = medium.sound_speed_m_s();
        let mach = velocity_m_s.map(|v| v / c);
        let speed = norm(mach);
        let displacement = core::array::from_fn(|i| {
            (Dd::from_f64(observer_m[i]) - Dd::from_f64(source_m[i])).to_f64()
        });
        let distance = norm(displacement);
        if !distance.is_finite() || distance < limits.minimum_distance_m {
            return Err(invalid("uniform-flow source/receiver separation"));
        }
        let along = dot(displacement.map(|v| v / distance), mach);
        let beta_squared = (1.0 - speed) * (1.0 + speed);
        let transformed = det::hypot(along, det::sqrt(beta_squared));
        // Rationalize the downstream root; subtracting two nearly equal
        // quantities loses the finite downstream delay as Mach approaches one.
        let factor = if along >= 0.0 {
            1.0 / (transformed + along)
        } else {
            (transformed - along) / beta_squared
        };
        let delay = (distance / c) * factor;
        if !delay.is_finite() || delay <= 0.0 {
            return Err(invalid("uniform-flow acoustic delay is unrepresentable"));
        }
        // Choose emission as the coordinate epoch. This avoids subtracting
        // U*t from large absolute positions on every receiver-time query.
        let receiver_in_rest_frame = core::array::from_fn(|i| {
            (Dd::from_f64(displacement[i]) - Dd::from_f64(velocity_m_s[i]) * Dd::from_f64(delay))
                .to_f64()
        });
        let rest_distance = norm(receiver_in_rest_frame);
        if !rest_distance.is_finite() || rest_distance <= 0.0 {
            return Err(invalid(
                "uniform-flow retarded separation is unrepresentable",
            ));
        }
        // Laboratory clearance was checked above. The rest-frame source event
        // has its own, generally shorter, nonzero distance; it is not a second
        // physical clearance criterion on the laboratory geometry.
        let pressure = MovingCompactObserver::new(
            medium,
            receiver_in_rest_frame,
            MovingCompactLimits {
                minimum_distance_m: rest_distance,
                ..limits
            },
        )?
        .prepare_emission(
            0.0,
            SourceKinematics {
                position_m: [0.0; 3],
                velocity_m_s: velocity_m_s.map(|v| -v),
                acceleration_m_s2: [0.0; 3],
            },
        )?;
        Ok(Self {
            receiver,
            source_m,
            pressure,
        })
    }

    /// Common receiver, retaining its fluid, laboratory velocity and limits.
    pub const fn receiver(&self) -> UniformFlowReceiver {
        self.receiver
    }
    /// Original fluid state; the flow never replaces its thermodynamic c or rho.
    pub const fn medium(&self) -> FluidAcousticState {
        self.receiver.medium
    }
    /// Uniform fluid velocity relative to the laboratory [m/s].
    pub const fn velocity_m_s(&self) -> [f64; 3] {
        self.receiver.velocity_m_s
    }
    /// Fixed physical source position [m].
    pub const fn source_position_m(&self) -> [f64; 3] {
        self.source_m
    }
    /// Fixed physical microphone position [m].
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.receiver.observer_m
    }
    /// Direction-dependent acoustic travel time [s], already used by phasors.
    pub const fn propagation_delay_s(&self) -> f64 {
        self.pressure.propagation_delay_s()
    }

    /// Prepare the actual emission time for one requested laboratory timestamp.
    ///
    /// The caller evaluates its original source history at `source_time_s()`.
    /// No extrapolation or resampling is performed. The represented timestamps
    /// must satisfy `|t-tau-T| <= maximum_time_residual_s`; unresolved arrivals
    /// at large epochs refuse instead of silently shifting the source clock.
    pub fn prepare_receiver(
        &self,
        receiver_time_s: f64,
        maximum_time_residual_s: f64,
    ) -> Result<PreparedUniformFlowSample, AeroacError> {
        if !receiver_time_s.is_finite()
            || !maximum_time_residual_s.is_finite()
            || maximum_time_residual_s < 0.0
        {
            return Err(invalid("uniform-flow receiver time and residual budget"));
        }
        let delay = self.propagation_delay_s();
        let source_time_s = (Dd::from_f64(receiver_time_s) - Dd::from_f64(delay)).to_f64();
        if !source_time_s.is_finite() || source_time_s >= receiver_time_s {
            return Err(invalid("uniform-flow emission time is unrepresentable"));
        }
        let error = residual(receiver_time_s, source_time_s, delay);
        if !error.is_finite() || error > maximum_time_residual_s {
            return Err(AeroacError::NotConverged {
                what: "uniform-flow retarded-time residual",
                residual: error,
            });
        }
        Ok(self.at_times(source_time_s, receiver_time_s, error))
    }

    /// Prepare an exact retained emission endpoint and its rounded arrival.
    /// The residual is reported, not waived: a history consumer must require
    /// the requested arrival to match and the residual to fit its time budget.
    pub fn prepare_emission(
        &self,
        source_time_s: f64,
    ) -> Result<PreparedUniformFlowSample, AeroacError> {
        let delay = self.propagation_delay_s();
        let receiver_time_s = (Dd::from_f64(source_time_s) + Dd::from_f64(delay)).to_f64();
        if !source_time_s.is_finite()
            || !receiver_time_s.is_finite()
            || receiver_time_s <= source_time_s
        {
            return Err(invalid(
                "uniform-flow emission or arrival time is unrepresentable",
            ));
        }
        let error = residual(receiver_time_s, source_time_s, delay);
        if !error.is_finite() {
            return Err(invalid("uniform-flow emission residual is unrepresentable"));
        }
        Ok(self.at_times(source_time_s, receiver_time_s, error))
    }

    fn at_times(
        &self,
        source_time_s: f64,
        receiver_time_s: f64,
        error: f64,
    ) -> PreparedUniformFlowSample {
        let mut pressure = self.pressure;
        pressure.sample.source_time_s = source_time_s;
        pressure.sample.receiver_time_s = receiver_time_s;
        pressure.sample.retarded_time_residual_s = error;
        PreparedUniformFlowSample {
            pressure,
            source_time_s,
            receiver_time_s,
            retarded_time_residual_s: error,
        }
    }

    /// Observe one original peak source phasor, using exp(-i omega t).
    ///
    /// This applies the same real pressure coefficients as transient queries,
    /// with exact harmonic time derivatives and the path phase exactly once.
    /// Wind changes direction-dependent amplitude and delay, but the frequency
    /// is unchanged for fixed laboratory endpoints. No plane-wave intensity or
    /// particle velocity is inferred from pressure in the moving medium.
    pub fn pressure_phasor(
        &self,
        frequency_hz: f64,
        source: CompactSourcePhasor,
    ) -> Result<C64, AeroacError> {
        let omega = core::f64::consts::TAU * frequency_hz;
        let phase = omega * self.propagation_delay_s();
        if !omega.is_finite() || omega <= 0.0 || !phase.is_finite() {
            return Err(invalid("uniform-flow acoustic frequency or phase"));
        }
        let q = source.volume_flow_m3_s;
        let f = source.force_on_fluid_n;
        let derivative = |v: C64| C64::new(omega * v.im, -omega * v.re);
        let values = [
            q,
            derivative(q),
            f[0],
            f[1],
            f[2],
            derivative(f[0]),
            derivative(f[1]),
            derivative(f[2]),
        ];
        let mut re = Dd::ZERO;
        let mut im = Dd::ZERO;
        for (coefficient, value) in self.pressure.coefficients.iter().zip(values) {
            if !finite(value) {
                return Err(nonfinite(
                    "uniform-flow harmonic source strength or derivative",
                ));
            }
            re = re + Dd::from_f64(*coefficient) * Dd::from_f64(value.re);
            im = im + Dd::from_f64(*coefficient) * Dd::from_f64(value.im);
        }
        let (sin, cos) = det::sin_cos(phase);
        let result = C64::new(
            (re * Dd::from_f64(cos) - im * Dd::from_f64(sin)).to_f64(),
            (re * Dd::from_f64(sin) + im * Dd::from_f64(cos)).to_f64(),
        );
        if !finite(result) {
            return Err(nonfinite("uniform-flow harmonic pressure"));
        }
        Ok(result)
    }
}

/// One admitted receiver/emission pair with the original pressure coefficients.
#[derive(Clone, Copy, Debug)]
pub struct PreparedUniformFlowSample {
    pressure: PreparedMovingCompactSample,
    source_time_s: f64,
    receiver_time_s: f64,
    retarded_time_residual_s: f64,
}

impl PreparedUniformFlowSample {
    /// Original Eulerian pressure coefficients with these physical timestamps.
    /// Galilean frames share time. The retained `source_time_rate()` is the
    /// fixed-position fluid-frame factor `1/kappa`, not the laboratory observed
    /// rate: fixed laboratory source and receiver have `d tau/d t = 1`.
    pub const fn pressure_in_fluid_frame(&self) -> PreparedMovingCompactSample {
        self.pressure
    }
    /// Original source time at which all strengths and derivatives are needed [s].
    pub const fn source_time_s(&self) -> f64 {
        self.source_time_s
    }
    /// Exact requested laboratory microphone time [s].
    pub const fn receiver_time_s(&self) -> f64 {
        self.receiver_time_s
    }
    /// Apply actual strengths at the admitted source time, without advancing it.
    pub fn evaluate(
        &self,
        strength: CompactSourceStrength,
    ) -> Result<UniformFlowSample, AeroacError> {
        let result = self.pressure.evaluate(strength)?;
        Ok(UniformFlowSample {
            source_time_s: self.source_time_s,
            receiver_time_s: self.receiver_time_s,
            propagation_delay_s: self.pressure.propagation_delay_s(),
            retarded_time_residual_s: self.retarded_time_residual_s,
            monopole_pressure_pa: result.monopole_pressure_pa,
            dipole_pressure_pa: result.dipole_pressure_pa,
            pressure_pa: result.pressure_pa,
        })
    }
}

/// Pressure in pascals on the original laboratory clock.
///
/// Since both laboratory endpoints are stationary, `d tau/d t = 1`.
/// The rest-frame stationary-receiver Doppler factor must not replace it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniformFlowSample {
    pub source_time_s: f64,
    pub receiver_time_s: f64,
    pub propagation_delay_s: f64,
    pub retarded_time_residual_s: f64,
    pub monopole_pressure_pa: f64,
    pub dipole_pressure_pa: f64,
    pub pressure_pa: f64,
}
