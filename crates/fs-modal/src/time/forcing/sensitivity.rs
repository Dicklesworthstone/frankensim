//! Physical-coefficient derivatives of the existing exact held-ramp flow.
//!
//! For k = omega^2 and c the physical viscous rate, the variational equations
//! are q_k'' + c q_k' + k q_k = -q and
//! q_c'' + c q_c' + k q_c = -q'. Initial sensitivities are zero. Differentiating
//! the entire flow moments and their exact composition solves these same
//! linear equations, without finite differences, pole separation or a special
//! critical-damping band. The step, source coordinates and force samples stay
//! fixed. This is a derivative of the linear dynamics, not a material law.

use super::{
    DampedTransition, LinearForcePullback, LinearForceTransition, ModalAcousticState,
    ModalAcousticTimeError, invalid, small_ramp_moments,
};

/// Endpoint derivatives with respect to the original physical coefficients.
///
/// `stiffness` means k = omega^2 [1/s^2], or K/m for a dimensional oscillator.
/// `damping_rate` means c [1/s], or C/m, held independent of k. Initial state,
/// endpoint forces, mass normalization and duration are fixed in every partial.
/// These fields are sensitivities, not displacements or velocities themselves.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearForceParameterSensitivity {
    /// dq_end/dk, in displacement units times s^2.
    pub displacement_by_stiffness: f64,
    /// dv_end/dk, in velocity units times s^2.
    pub velocity_by_stiffness: f64,
    /// dq_end/dc, in displacement units times s.
    pub displacement_by_damping_rate: f64,
    /// dv_end/dc, in velocity units times s.
    pub velocity_by_damping_rate: f64,
}

/// Original endpoint motion and its physical-coefficient sensitivities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearForceSensitivityFrame {
    /// The original [`LinearForceTransition::advance`] result, unchanged.
    pub state: ModalAcousticState,
    /// Fixed-input derivatives of that continuous oscillator flow.
    pub parameters: LinearForceParameterSensitivity,
}

/// Cotangents of all state, force and physical-coefficient inputs to one step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearForceParameterPullback {
    /// Initial-state and force-sample cotangents from the original transpose.
    pub inputs: LinearForcePullback,
    /// dJ/dk for k = omega^2 [1/s^2], holding c fixed.
    pub stiffness: f64,
    /// dJ/dc for the physical damping rate c [1/s], holding k fixed.
    pub damping_rate: f64,
}

/// Prepared exact physical damping/stiffness sensitivities for a fixed step.
///
/// Preparation differentiates the same bounded small-interval entire series
/// and algebraic doubling identities as the original flow. It does not
/// differentiate an iterative equation solver or perturb complete trajectories.
/// Runtime advancement and pullback allocate nothing, call no transcendental
/// functions, and retain the original forward step's arithmetic.
///
/// k and c remain independent, including k = 0 and critical damping. For an
/// omega parameter at fixed c, apply dJ/domega = 2*omega*dJ/dk. For a ratio
/// parameterization c = 2*zeta*omega, the chain rule instead gives
/// dJ/dzeta = 2*omega*dJ/dc and
/// dJ/domega = 2*omega*dJ/dk + 2*zeta*dJ/dc. Source/material and changing-basis
/// derivatives remain the caller's responsibility. Work/heat integrals,
/// duration and nonsmooth observations are not differentiated by this owner.
#[derive(Clone, Copy, Debug)]
pub struct PreparedLinearForceSensitivity {
    motion: LinearForceTransition,
    stiffness: MomentDerivative,
    damping: MomentDerivative,
}

#[derive(Clone, Copy, Debug, Default)]
struct MomentDerivative {
    b: f64,
    g: f64,
    ramp: f64,
    position: f64,
    velocity: f64,
}

#[derive(Clone, Copy)]
enum Parameter {
    Stiffness,
    Damping,
}

impl Parameter {
    fn seeds(self) -> (f64, f64, f64) {
        // Local variables are A = (omega*h)^2 and d = c*h. Doubling h
        // multiplies them by four/two, so derivatives transform inversely.
        match self {
            Self::Stiffness => (1.0, 0.0, 0.25),
            Self::Damping => (0.0, 1.0, 0.5),
        }
    }
}

impl LinearForceTransition {
    /// Prepare analytic stiffness and damping-rate derivatives for this step.
    ///
    /// The original forward map is unchanged. The extra scalar preparation is
    /// opt-in and reusable across states, loads and reverse trajectory sweeps.
    /// A ratio-based original step is differentiated at fixed physical c, not
    /// fixed zeta; see [`PreparedLinearForceSensitivity`] for the chain rules.
    ///
    /// # Errors
    /// Unrepresentable dimensionless derivative coefficients. Runtime methods
    /// separately refuse unrepresentable input-dependent physical derivatives.
    pub fn prepare_sensitivities(
        &self,
    ) -> Result<PreparedLinearForceSensitivity, ModalAcousticTimeError> {
        let mut angle = self.transition.angle;
        let mut drag = self.transition.drag;
        let mut doublings = 0;
        while angle.max(drag) > 0.25 {
            angle *= 0.5;
            drag *= 0.5;
            doublings += 1;
        }
        let mut flow = DampedTransition::small(angle, drag);
        let mut ramp = small_ramp_moments(angle, drag).0;
        let mut stiffness = MomentDerivative::small(flow, Parameter::Stiffness);
        let mut damping = MomentDerivative::small(flow, Parameter::Damping);
        for _ in 0..doublings {
            let before = flow;
            flow.double();
            stiffness = stiffness.double(before, flow, ramp, Parameter::Stiffness);
            damping = damping.double(before, flow, ramp, Parameter::Damping);
            ramp = ((1.0 + before.position_retention) * ramp + (1.0 + before.b) * before.g) * 0.125;
        }
        if !stiffness.is_finite() || !damping.is_finite() {
            return Err(invalid(
                "linear-force coefficient sensitivities must be representable",
            ));
        }
        Ok(PreparedLinearForceSensitivity {
            motion: *self,
            stiffness,
            damping,
        })
    }
}

impl PreparedLinearForceSensitivity {
    /// Advance the original motion and evaluate both coefficient sensitivities.
    ///
    /// The initial state is independent of the parameters for this local
    /// partial. A trajectory sensitivity must also propagate its incoming
    /// state derivative; [`Self::pullback`] supplies the reverse-mode form.
    ///
    /// # Errors
    /// Every refusal from [`LinearForceTransition::advance`], or a nonfinite
    /// physical sensitivity. No state is mutated or partial result returned.
    pub fn advance(
        &self,
        state: ModalAcousticState,
        force_start: f64,
        force_end: f64,
    ) -> Result<LinearForceSensitivityFrame, ModalAcousticTimeError> {
        let (next, residual, change) = self.motion.motion(state, force_start, force_end)?;
        let h = self.motion.duration_s;
        let flow = self.motion.transition;
        let q = state.displacement_m_sqrt_kg;
        let v = state.velocity_m_sqrt_kg_per_s;
        let k = self.stiffness;
        let c = self.damping;
        let qk = if flow.position_decrement.abs() < 0.5 {
            // Differentiate the original imbalance, not two large equilibrium
            // positions: d[h*(F-k*q)]/dk = -h*q. Group h factors so h^2
            // need not overflow before a small derivative rescales it.
            h * (h * (h * (k.b * v + k.g * residual + k.ramp * change) - flow.g * q))
        } else {
            // Retain the independently differentiated homogeneous tail when
            // forming it as one minus a decrement would cancel it away.
            let force = h * force_start;
            h * (h * (k.position * q + h * (k.b * v + k.g * force + k.ramp * change)))
        };
        // The damping force moment retains the differentiated homogeneous
        // tail. Its imbalance form also makes every exactly balanced static
        // state independent of c, without subtracting two parameter responses.
        let qc = h * (h * (c.b * v + c.g * residual + c.ramp * change));
        let vk = h * (h * (k.velocity * v + k.b * residual + k.g * change) - flow.b * q);
        let vc = h * (c.velocity * v + c.b * residual + c.g * change);
        if ![qk, qc, vk, vc].into_iter().all(f64::is_finite) {
            return Err(invalid(
                "linear-force physical sensitivities must be representable",
            ));
        }
        Ok(LinearForceSensitivityFrame {
            state: next,
            parameters: LinearForceParameterSensitivity {
                displacement_by_stiffness: qk,
                velocity_by_stiffness: vk,
                displacement_by_damping_rate: qc,
                velocity_by_damping_rate: vc,
            },
        })
    }

    /// Pull endpoint cotangents through the original state/force/coefficient map.
    ///
    /// Add both coefficient cotangents over all steps sharing the same k/c,
    /// propagate both initial-state cotangents, and sum the two contributions
    /// at each shared force sample. Dependence of the initial state or forcing
    /// on the material parameters requires its own chain rule.
    ///
    /// # Errors
    /// Invalid state, force or cotangent inputs, or unrepresentable motion,
    /// sensitivity or pullback arithmetic. The prepared flow is never mutated.
    pub fn pullback(
        &self,
        state: ModalAcousticState,
        force_start: f64,
        force_end: f64,
        endpoint_displacement_cotangent: f64,
        endpoint_velocity_cotangent: f64,
    ) -> Result<LinearForceParameterPullback, ModalAcousticTimeError> {
        let inputs = self
            .motion
            .pullback(endpoint_displacement_cotangent, endpoint_velocity_cotangent)?;
        let parameter = self.advance(state, force_start, force_end)?.parameters;
        let stiffness = endpoint_displacement_cotangent.mul_add(
            parameter.displacement_by_stiffness,
            endpoint_velocity_cotangent * parameter.velocity_by_stiffness,
        );
        let damping_rate = endpoint_displacement_cotangent.mul_add(
            parameter.displacement_by_damping_rate,
            endpoint_velocity_cotangent * parameter.velocity_by_damping_rate,
        );
        if ![stiffness, damping_rate].into_iter().all(f64::is_finite) {
            return Err(invalid(
                "linear-force coefficient cotangents must be representable",
            ));
        }
        Ok(LinearForceParameterPullback {
            inputs,
            stiffness,
            damping_rate,
        })
    }
}

impl MomentDerivative {
    fn small(flow: DampedTransition, parameter: Parameter) -> Self {
        let (da, dd, _) = parameter.seeds();
        let a = flow.angle * flow.angle;
        let d = flow.drag;
        let mut previous = 0.0;
        let mut power = 1.0;
        let mut previous_derivative = 0.0;
        let mut power_derivative = 0.0;
        let mut reciprocal_factorial = 1.0;
        let mut result = Self::default();
        for j in 1..=18 {
            let next = -d * power - a * previous;
            let derivative =
                -dd * power - d * power_derivative - da * previous - a * previous_derivative;
            previous = power;
            power = next;
            previous_derivative = power_derivative;
            power_derivative = derivative;
            reciprocal_factorial /= f64::from(j + 1);
            let term = power_derivative * reciprocal_factorial;
            result.b += term;
            result.g += term / f64::from(j + 2);
            result.ramp += term / (f64::from(j + 2) * f64::from(j + 3));
        }
        result.position = -da * flow.g - a * result.g;
        result.velocity = result.position - dd * flow.b - d * result.b;
        result
    }

    fn double(
        self,
        before: DampedTransition,
        after: DampedTransition,
        ramp: f64,
        parameter: Parameter,
    ) -> Self {
        let (da, dd, scale) = parameter.seeds();
        let cross =
            da * before.b * before.b + 2.0 * (before.angle * before.b) * (before.angle * self.b);
        let mut next = Self {
            b: 0.5
                * scale
                * (self.b * (before.position_retention + before.velocity_retention)
                    + before.b * (self.position + self.velocity)),
            g: 0.25
                * scale
                * (self.g * (1.0 + before.position_retention)
                    + before.g * self.position
                    + 2.0 * before.b * self.b),
            ramp: 0.125
                * scale
                * ((1.0 + before.position_retention) * self.ramp
                    + ramp * self.position
                    + (1.0 + before.b) * self.g
                    + before.g * self.b),
            position: scale * (2.0 * before.position_retention * self.position - cross),
            velocity: scale * (2.0 * before.velocity_retention * self.velocity - cross),
        };
        let position = -da * after.g - after.angle * (after.angle * next.g);
        if after.position_decrement.abs() < 0.5 {
            next.position = position;
        } else if matches!(parameter, Parameter::Damping) {
            // P = 1 - A*g gives g_d = -P_d/A. Once the position change
            // is large, A is safely nonzero. Preserve the decaying forced
            // response derivative instead of cancelling it inside g's
            // composition at equilibrium. Scale the quotient in two stages.
            next.g = -(next.position / after.angle) / after.angle;
        }
        if after.velocity_decrement.abs() < 0.5 {
            next.velocity = next.position - dd * after.b - after.drag * next.b;
        }
        next
    }

    fn is_finite(self) -> bool {
        [self.b, self.g, self.ramp, self.position, self.velocity]
            .into_iter()
            .all(f64::is_finite)
    }
}
