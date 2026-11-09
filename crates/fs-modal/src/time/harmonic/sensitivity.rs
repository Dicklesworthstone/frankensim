//! Physical coefficient tangents of the original continuous-tone convolution.
//!
//! At fixed forcing clock, q_k'' + c q_k' + k q_k = -q and
//! q_c'' + c q_c' + k q_c = -q'. Differentiating the entire small-step
//! action and its doubling solves these variational equations through startup
//! and resonance, without finite differences or a steady-state approximation.

use super::{
    DampedTransition, HarmonicForceTransition, ModalAcousticState, ModalAcousticTimeError,
    double_force, invalid, small_force,
};
use crate::time::forcing::{
    sensitivity::{MomentDerivative, Parameter},
    small_ramp_moments,
};
use fs_math::det;

mod clock;
pub use clock::HarmonicForcePhasePullback;

/// Endpoint partials in the original mass-normalized physical coordinates.
/// Each partial holds the initial state, duration and other physical inputs
/// fixed. Initial phase means a common additional phase rotation of the two
/// prescribed force quadratures at the interval start.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicForceParameterSensitivity {
    /// dq_end/dk for k = omega^2 [1/s^2], holding c fixed.
    pub displacement_by_stiffness: f64,
    /// dv_end/dk in velocity units times s^2.
    pub velocity_by_stiffness: f64,
    /// dq_end/dc for physical viscous rate c [1/s], holding k fixed.
    pub displacement_by_damping_rate: f64,
    /// dv_end/dc in velocity units times s.
    pub velocity_by_damping_rate: f64,
    /// dq_end/dOmega, varying the FORCING angular frequency [rad/s] at fixed
    /// initial source phase. This is independent of natural-frequency/k changes.
    pub displacement_by_forcing_angular_frequency: f64,
    /// dv_end/dOmega in velocity units times seconds.
    pub velocity_by_forcing_angular_frequency: f64,
    /// dq_end/dphi for F = Fc*cos(Omega*t+phi) + Fs*sin(Omega*t+phi), at phi=0.
    pub displacement_by_initial_phase: f64,
    /// dv_end/dphi for the same common initial phase [rad].
    pub velocity_by_initial_phase: f64,
}

/// Original motion plus its independently propagated physical partials.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicForceSensitivityFrame {
    /// The original [`HarmonicForceTransition::advance`] result, unchanged.
    pub state: ModalAcousticState,
    /// Local physical stiffness, viscous-rate and source clock/phase derivatives.
    pub parameters: HarmonicForceParameterSensitivity,
}

/// Scalar-objective derivatives with respect to this step's physical inputs.
/// Cosine and sine forces refer to the SAME interval-start phase as the motion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicForcePullback {
    /// dJ/dq_start.
    pub initial_displacement: f64,
    /// dJ/dv_start.
    pub initial_velocity: f64,
    /// dJ/dFc in the interval-start cosine basis.
    pub cosine_force: f64,
    /// dJ/dFs in the interval-start sine basis.
    pub sine_force: f64,
    /// dJ/dk for k = omega^2, at fixed c.
    pub stiffness: f64,
    /// dJ/dc for the physical viscous rate, at fixed k.
    pub damping_rate: f64,
    /// dJ/dOmega for the FORCING angular frequency, at fixed initial phase.
    pub forcing_angular_frequency: f64,
    /// dJ/dphi for a common additional interval-start source phase [rad].
    pub initial_phase: f64,
}

/// Analytic harmonic-transient derivatives prepared outside the hot loop.
///
/// The original entire convolution is differentiated, using the existing
/// homogeneous-flow derivative owner. No detuning or pole separation is
/// divided out: free, critical, overdamped and exactly resonant modes retain
/// one formula. Runtime advance/pullback allocate nothing and evaluate no
/// transcendental functions. Forward motion keeps the original cached bits.
///
/// Partials use independent k = omega^2 and c = C/m, including their one-sided
/// limits at zero. At fixed c, dJ/domega = 2*omega*dJ/dk. For c=2*zeta*omega,
/// dJ/domega also includes 2*zeta*dJ/dc, and dJ/dzeta=2*omega*dJ/dc.
/// Forcing-frequency derivatives instead integrate the exact affine envelope
/// t*(Fs*cos(Omega*t)-Fc*sin(Omega*t)); phase derivatives rotate the two force
/// quadratures. Duration, material/basis changes and work/heat are not
/// differentiated. A reverse trajectory propagates state cotangents and sums
/// shared parameter contributions. When it advances the source phase too,
/// [`HarmonicForceTransition::pullback_force_phase`] supplies that clock chain
/// rule; a local frequency partial alone holds each interval's input phase fixed.
#[derive(Clone, Copy, Debug)]
pub struct PreparedHarmonicForceSensitivity {
    motion: HarmonicForceTransition,
    stiffness: MomentDerivative,
    damping: MomentDerivative,
    force_stiffness: [f64; 4],
    force_damping: [f64; 4],
    force_envelope: [f64; 4],
}

impl HarmonicForceTransition {
    /// Cold-prepare stiffness, physical damping-rate and prescribed source
    /// clock/phase derivatives. The original forward cache and state are read-only.
    ///
    /// # Errors
    /// Unrepresentable dimensionless derivative coefficients. Runtime queries
    /// separately refuse input-dependent physical overflow without a partial result.
    pub fn prepare_sensitivities(
        &self,
    ) -> Result<PreparedHarmonicForceSensitivity, ModalAcousticTimeError> {
        let (mut angle, mut drag, mut phase) = (
            self.transition.angle,
            self.transition.drag,
            self.forcing_phase_rad,
        );
        let mut doublings = 0;
        while angle.max(drag).max(phase) > 0.25 {
            angle *= 0.5;
            drag *= 0.5;
            phase *= 0.5;
            doublings += 1;
        }
        let mut flow = DampedTransition::small(angle, drag);
        let mut ramp = small_ramp_moments(angle, drag).0;
        let mut stiffness = MomentDerivative::small(flow, Parameter::Stiffness);
        let mut damping = MomentDerivative::small(flow, Parameter::Damping);
        let remainder = self.cosine_remainder;
        let mut force = small_force(angle, drag, phase, remainder);
        let mut force_stiffness =
            small_derivative(angle, drag, phase, remainder, Parameter::Stiffness);
        let mut force_damping = small_derivative(angle, drag, phase, remainder, Parameter::Damping);
        for _ in 0..doublings {
            force_stiffness = double_derivative(
                force_stiffness,
                force,
                flow,
                stiffness,
                phase,
                remainder,
                Parameter::Stiffness,
            );
            force_damping = double_derivative(
                force_damping,
                force,
                flow,
                damping,
                phase,
                remainder,
                Parameter::Damping,
            );
            force = double_force(force, flow, phase, remainder);
            let before = flow;
            flow.double();
            stiffness = stiffness.double(before, flow, ramp, Parameter::Stiffness);
            damping = damping.double(before, flow, ramp, Parameter::Damping);
            ramp = ((1.0 + before.position_retention) * ramp + (1.0 + before.b) * before.g) * 0.125;
            phase *= 2.0;
        }
        if !stiffness.is_finite()
            || !damping.is_finite()
            || force_stiffness
                .into_iter()
                .chain(force_damping)
                .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "harmonic coefficient sensitivities must be representable",
            ));
        }
        Ok(PreparedHarmonicForceSensitivity {
            motion: *self,
            stiffness,
            damping,
            force_stiffness,
            force_damping,
            force_envelope: super::affine::delta_gains_for_phase(self, self.forcing_phase_rad)?,
        })
    }
}

impl PreparedHarmonicForceSensitivity {
    /// Advance the original motion and evaluate local fixed-input partials.
    /// Incoming trajectory state sensitivities are not assumed zero globally:
    /// propagate them through the homogeneous map, or use [`Self::pullback`].
    ///
    /// # Errors
    /// Original motion refusal or unrepresentable physical derivatives.
    pub fn advance(
        &self,
        state: ModalAcousticState,
        cosine_force: f64,
        sine_force: f64,
    ) -> Result<HarmonicForceSensitivityFrame, ModalAcousticTimeError> {
        let next = self.motion.advance(state, cosine_force, sine_force)?;
        let h = self.motion.duration_s;
        let flow = self.motion.transition;
        let q = state.displacement_m_sqrt_kg;
        let v = state.velocity_m_sqrt_kg_per_s;
        let fc = h * cosine_force;
        let fs = h * sine_force;
        let held = if self.motion.cosine_remainder {
            fc
        } else {
            0.0
        };
        let restoring = (self.motion.angular_frequency_rad_s * q) * flow.angle;
        let residual = held - restoring;
        let k = self.stiffness;
        let c = self.damping;
        let [kqc, kqs, kvc, kvs] = self.force_stiffness;
        let [cqc, cqs, cvc, cvs] = self.force_damping;
        let force_qk = kqc * fc + kqs * fs;
        // Retain the original force imbalance around static equilibrium and
        // independently computed homogeneous tails after a long decay.
        let qk = if flow.position_decrement.abs() < 0.5 {
            h * (h * (h * (k.b * v + k.g * residual + force_qk) - flow.g * q))
        } else {
            h * (h * (k.position * q + h * (k.b * v + k.g * held + force_qk)))
        };
        let qc = h * (h * (c.b * v + c.g * residual + cqc * fc + cqs * fs));
        let vk = h * (h * (k.velocity * v + k.b * residual + kvc * fc + kvs * fs) - flow.b * q);
        let vc = h * (c.velocity * v + c.b * residual + cvc * fc + cvs * fs);
        if ![qk, vk, qc, vc].into_iter().all(f64::is_finite) {
            return Err(invalid(
                "harmonic physical sensitivities must be representable",
            ));
        }
        let [qw, vw, qp, vp] = self.clock_partials(cosine_force, sine_force)?;
        Ok(HarmonicForceSensitivityFrame {
            state: next,
            parameters: HarmonicForceParameterSensitivity {
                displacement_by_stiffness: qk,
                velocity_by_stiffness: vk,
                displacement_by_damping_rate: qc,
                velocity_by_damping_rate: vc,
                displacement_by_forcing_angular_frequency: qw,
                velocity_by_forcing_angular_frequency: vw,
                displacement_by_initial_phase: qp,
                velocity_by_initial_phase: vp,
            },
        })
    }

    /// Pull endpoint-objective cotangents through state, both force quadratures,
    /// stiffness, physical viscous rate and source clock/phase on the same
    /// continuous interval. Initial phase is independent of forcing frequency.
    ///
    /// # Errors
    /// Invalid state, force or cotangents, or unrepresentable motion/derivatives.
    pub fn pullback(
        &self,
        state: ModalAcousticState,
        cosine_force: f64,
        sine_force: f64,
        endpoint_displacement_cotangent: f64,
        endpoint_velocity_cotangent: f64,
    ) -> Result<HarmonicForcePullback, ModalAcousticTimeError> {
        let (q_bar, v_bar) = (endpoint_displacement_cotangent, endpoint_velocity_cotangent);
        if ![q_bar, v_bar].into_iter().all(f64::is_finite) {
            return Err(invalid("harmonic endpoint cotangents must be finite"));
        }
        let parameter = self.advance(state, cosine_force, sine_force)?.parameters;
        let h = self.motion.duration_s;
        let flow = self.motion.transition;
        let restoring = self.motion.angular_frequency_rad_s * ((flow.angle * flow.b) * v_bar);
        let initial_displacement = if flow.position_decrement.abs() < 0.5 {
            q_bar - (flow.position_decrement * q_bar + restoring)
        } else {
            flow.position_retention * q_bar - restoring
        };
        let scaled_q_bar = h * q_bar;
        let initial_velocity = if flow.velocity_decrement.abs() < 0.5 {
            v_bar + (flow.b * scaled_q_bar - flow.velocity_decrement * v_bar)
        } else {
            flow.b * scaled_q_bar + flow.velocity_retention * v_bar
        };
        let [mut qc, qs, mut vc, vs] = self.motion.force;
        if self.motion.cosine_remainder {
            qc += flow.g;
            vc += flow.b;
        }
        let cosine_force = h * (qc * scaled_q_bar + vc * v_bar);
        let sine_force = h * (qs * scaled_q_bar + vs * v_bar);
        let stiffness = q_bar.mul_add(
            parameter.displacement_by_stiffness,
            v_bar * parameter.velocity_by_stiffness,
        );
        let damping_rate = q_bar.mul_add(
            parameter.displacement_by_damping_rate,
            v_bar * parameter.velocity_by_damping_rate,
        );
        let forcing_angular_frequency = q_bar.mul_add(
            parameter.displacement_by_forcing_angular_frequency,
            v_bar * parameter.velocity_by_forcing_angular_frequency,
        );
        let initial_phase = q_bar.mul_add(
            parameter.displacement_by_initial_phase,
            v_bar * parameter.velocity_by_initial_phase,
        );
        if ![
            initial_displacement,
            initial_velocity,
            cosine_force,
            sine_force,
            stiffness,
            damping_rate,
            forcing_angular_frequency,
            initial_phase,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(invalid("harmonic input cotangents must be representable"));
        }
        Ok(HarmonicForcePullback {
            initial_displacement,
            initial_velocity,
            cosine_force,
            sine_force,
            stiffness,
            damping_rate,
            forcing_angular_frequency,
            initial_phase,
        })
    }
}

fn small_derivative(
    angle: f64,
    drag: f64,
    phase: f64,
    remainder: bool,
    parameter: Parameter,
) -> [f64; 4] {
    let (da, dd, _) = parameter.seeds();
    let mut gains = [0.0; 4];
    for column in 0..2 {
        let mut term = [0.0; 4];
        term[2 + column] = 1.0;
        let mut derivative = [0.0; 2];
        for order in 1..=24 {
            let inverse = 1.0 / f64::from(order);
            derivative = [
                derivative[1] * inverse,
                (-da * term[0]
                    - angle * angle * derivative[0]
                    - dd * term[1]
                    - drag * derivative[1])
                    * inverse,
            ];
            let source = if remainder && column == 0 && order == 1 {
                0.0
            } else {
                term[2]
            };
            term = [
                term[1] * inverse,
                (-angle * angle * term[0] - drag * term[1] + source) * inverse,
                phase * term[3] * inverse,
                -phase * term[2] * inverse,
            ];
            gains[column] += derivative[0];
            gains[2 + column] += derivative[1];
        }
    }
    gains
}

fn double_derivative(
    derivative: [f64; 4],
    force: [f64; 4],
    flow: DampedTransition,
    tangent: MomentDerivative,
    phase: f64,
    remainder: bool,
    parameter: Parameter,
) -> [f64; 4] {
    let (da, _, scale) = parameter.seeds();
    let (sine, cosine) = det::sin_cos(phase);
    let [qc, qs, vc, vs] = force;
    let [dqc, dqs, dvc, dvs] = derivative;
    let restoring = -(flow.angle * flow.b) * flow.angle;
    let dr = -da * flow.b - flow.angle * (flow.angle * tangent.b);
    let (constant_cos, constant_sin) = if remainder {
        let half_sine = det::sin(0.5 * phase);
        (-2.0 * half_sine * half_sine, sine)
    } else {
        (0.0, 0.0)
    };
    [
        0.25 * scale
            * (tangent.position * qc
                + flow.position_retention * dqc
                + tangent.b * vc
                + flow.b * dvc
                + dqc * cosine
                - dqs * sine
                + tangent.g * constant_cos),
        0.25 * scale
            * (tangent.position * qs
                + flow.position_retention * dqs
                + tangent.b * vs
                + flow.b * dvs
                + dqc * sine
                + dqs * cosine
                + tangent.g * constant_sin),
        0.5 * scale
            * (dr * qc
                + restoring * dqc
                + tangent.velocity * vc
                + flow.velocity_retention * dvc
                + dvc * cosine
                - dvs * sine
                + tangent.b * constant_cos),
        0.5 * scale
            * (dr * qs
                + restoring * dqs
                + tangent.velocity * vs
                + flow.velocity_retention * dvs
                + dvc * sine
                + dvs * cosine
                + tangent.b * constant_sin),
    ]
}
