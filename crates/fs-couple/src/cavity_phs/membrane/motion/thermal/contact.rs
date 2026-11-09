//! Reciprocal finite-body heat contact on the original nonlinear membrane step.
//!
//! The native Gibbs pressure drives membrane momentum and volume work. The
//! original body enthalpy chart supplies its entropy-conjugate temperature.
//! A single signed heat Q enters the first gas and leaves the finite body.

use super::{
    GasMembraneDiscreteGradientStep, GasMembraneMotionError, GasMembraneThermalControls,
    GasMembraneThermalMotion,
};
use crate::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve,
    feedback::{self, ScalarError},
};
use fs_material::phase::{
    EquilibriumPhaseState, PhaseStateError, UniformEnthalpyStep, UniformEnthalpyStepInput,
};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::parcel::discrete_gradient::GasDiscreteGradientStep;
use fs_material::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;
use fs_math::dd::Dd;
use fs_time::galpha::SecondOrderState;

/// The actual original source state used to sample thermal contact.
/// Its EOS pressure is distinct from the discrete mechanical work pressure.
#[derive(Clone, Debug, PartialEq)]
pub enum MembraneGasSample {
    Ideal(ResolvedThermochemicalGasState),
    PengRobinson(ResolvedThermochemicalRealGasState),
}
impl MembraneGasSample {
    /// Entropy-conjugate source temperature [K].
    pub fn temperature_k(&self) -> f64 {
        match self {
            Self::Ideal(state) => state.state().temperature,
            Self::PengRobinson(state) => state.temperature_k(),
        }
    }
    /// Actual EOS pressure of the source sample [Pa].
    pub fn pressure_pa(&self) -> f64 {
        match self {
            Self::Ideal(state) => state.state().pressure,
            Self::PengRobinson(state) => state.pressure_pa(),
        }
    }
}

/// Original physical refusals and bounded heat-contact failures.
#[derive(Debug)]
pub enum GasMembraneContactError<E> {
    Invalid(&'static str),
    Motion(GasMembraneMotionError),
    Thermal(PhaseStateError),
    Interface(E),
    Cancelled,
    OutsideBracket {
        lower_residual_j: f64,
        upper_residual_j: f64,
    },
    NotConverged {
        iterations: usize,
        residual_j: f64,
    },
}
impl<E> From<GasMembraneMotionError> for GasMembraneContactError<E> {
    fn from(error: GasMembraneMotionError) -> Self {
        match error {
            GasMembraneMotionError::Cancelled => Self::Cancelled,
            other => Self::Motion(other),
        }
    }
}
impl<E> From<PhaseStateError> for GasMembraneContactError<E> {
    fn from(error: PhaseStateError) -> Self {
        match error {
            PhaseStateError::EntropyEvaluationCancelled => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}
impl<E: core::fmt::Display> core::fmt::Display for GasMembraneContactError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "gas/membrane/body contact: {reason}"),
            Self::Motion(error) => write!(f, "gas/membrane/body contact: {error}"),
            Self::Thermal(error) => write!(f, "gas/membrane/body material: {error}"),
            Self::Interface(error) => write!(f, "gas/membrane/body interface: {error}"),
            Self::Cancelled => write!(f, "gas/membrane/body contact cancelled"),
            Self::OutsideBracket {
                lower_residual_j,
                upper_residual_j,
            } => write!(
                f,
                "membrane heat contact not bracketed: {lower_residual_j}, {upper_residual_j} J"
            ),
            Self::NotConverged {
                iterations,
                residual_j,
            } => write!(
                f,
                "membrane heat contact did not converge in {iterations} trials: {residual_j} J"
            ),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for GasMembraneContactError<E> {}

/// Independently evaluated mechanics, finite-body energy and thermal entropy.
#[derive(Clone, Debug)]
pub struct GasMembraneContactStep {
    /// Native trajectory, Gibbs gas receipts and separate pressure-work defect.
    pub motion: GasMembraneDiscreteGradientStep,
    /// Original source state retained from the accepted mechanical trial.
    pub gas_thermodynamic_mean: MembraneGasSample,
    /// Original body state at midpoint enthalpy, retaining its actual phase.
    pub thermal_midpoint: EquilibriumPhaseState,
    /// Original body Delta h / Delta s [K], including latent intervals.
    pub thermal_mean_temperature_k: f64,
    pub conductance_w_k: f64,
    /// Signed heat into the first gas and out of the body [J].
    pub exchanged_heat_j: f64,
    /// Independently sampled dt G (theta_body - theta_gas) [J].
    pub evaluated_heat_j: f64,
    pub source_residual_j: f64,
    pub thermal_energy_change_j: f64,
    pub thermal_energy_residual_j: f64,
    /// All gas, body, membrane and kinetic storage plus exterior work minus
    /// dead-load work and prescribed heat. Reciprocal Q cancels internally.
    pub total_energy_residual_j: f64,
    pub thermal_entropy_change_j_k: f64,
    /// Original gas(es) plus body entropy; structural integration error remains
    /// in motion.integration_defect_j and is never converted into entropy.
    pub total_entropy_change_j_k: f64,
    /// Prescribed gas/body/opposing-gas heat divided by each conjugate
    /// temperature. The external sources' own entropy remains caller-owned.
    pub prescribed_source_entropy_j_k: f64,
    /// Evaluated dt G (theta_body - theta_gas)^2 / (theta_body theta_gas).
    pub contact_entropy_production_j_k: f64,
    pub entropy_residual_j_k: f64,
    /// Allowance derived from original source, body and heat-law solve budgets.
    pub entropy_tolerance_j_k: f64,
    /// Complete native mechanical/contact proposals, including discarded ones.
    pub exchange_trials: usize,
    /// Native operator evaluations summed over all contact proposals.
    pub motion_evaluations: usize,
    /// Native gas/geometry trials summed over all contact proposals.
    pub gas_trials: usize,
}

/// Joint motion and finite-body proposal. The original accepted inputs are
/// unchanged; publish both returned states in the enclosing physical owner.
#[derive(Clone)]
pub struct GasMembraneContactCandidate<'m, 'p, S> {
    motion: GasMembraneThermalMotion<'m, 'p>,
    thermal: UniformEnthalpyStep<S>,
    step: GasMembraneContactStep,
}
impl<'m, 'p, S> GasMembraneContactCandidate<'m, 'p, S> {
    pub fn motion(&self) -> &GasMembraneThermalMotion<'m, 'p> {
        &self.motion
    }
    pub fn thermal(&self) -> &UniformEnthalpyStep<S> {
        &self.thermal
    }
    pub fn step(&self) -> &GasMembraneContactStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        GasMembraneThermalMotion<'m, 'p>,
        UniformEnthalpyStep<S>,
        GasMembraneContactStep,
    ) {
        (self.motion, self.thermal, self.step)
    }
}

fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasMembraneContactError<E>> {
    if cancelled() {
        Err(GasMembraneContactError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

impl<'m, 'p> GasMembraneThermalMotion<'m, 'p> {
    /// Solve finite-body thermal contact and nonlinear membrane pressure motion
    /// over the same physical interval, using the original discrete Gibbs law.
    ///
    /// The body exchanges with the first gas; a sealed opposing gas retains its
    /// separate prescribed heat. All three prescribed heat inputs exclude Q.
    /// The contact law is Q = dt G (theta_body - theta_gas), with each theta
    /// supplied by its original entropy/energy chart. The callback receives the
    /// native gas thermodynamic sample and body enthalpy midpoint, preserving
    /// phase-dependent transport. It must be repeatable and defer side effects.
    ///
    /// Source/geometry admission, native pressure-work and mechanical-energy
    /// gates remain active in every root trial. One scratch trajectory is reused;
    /// accepted history is copied only for the final joint candidate. Existing
    /// midpoint and prescribed-heat Gibbs methods retain their public behavior.
    ///
    /// The finite body has fixed mass/geometry and the original negligible-pV-work
    /// enthalpy storage model. It is an explicit thermal inventory, not inferred
    /// membrane mass or thermoelastic expansion. Uniform-body validity, contact
    /// geometry and conductance fidelity belong to its source. Positive contact
    /// production does not certify temporal accuracy; refine the physical step.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_discrete_gradient_candidate<S, E>(
        &self,
        duration_s: f64,
        dead_load_n: &[[f64; 3]],
        heat_j: f64,
        exterior_heat_j: Option<f64>,
        body: UniformEnthalpyStepInput<'_>,
        controls: GasMembraneThermalControls,
        exchange: HeatExchangeSolve,
        mut resolve: impl FnMut(
            &MembraneGasSample,
            EquilibriumPhaseState,
        ) -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasMembraneContactCandidate<'m, 'p, S>, GasMembraneContactError<E>> {
        poll(&mut cancelled)?;
        feedback::admit::<()>(exchange).map_err(|_| {
            GasMembraneContactError::Invalid("invalid finite-body heat bracket or work limits")
        })?;
        controls.validate()?;
        if !positive(body.mass_kg)
            || !positive(body.volume_m3)
            || !positive(body.surface_area_m2)
            || !positive(body.duration_s)
            || body.duration_s != duration_s
            || !body.internal_heat_j.is_finite()
            || !heat_j.is_finite()
            || exterior_heat_j.is_some_and(|q| !q.is_finite())
        {
            return Err(GasMembraneContactError::Invalid(
                "positive finite body geometry, matching duration and finite prescribed heat required",
            ));
        }
        let h0 = body.initial.specific_enthalpy_j_kg();
        if body.curve.state_at_specific_enthalpy(h0)? != body.initial {
            return Err(GasMembraneContactError::Invalid(
                "initial body state belongs to a different enthalpy chart",
            ));
        }
        // Reuse trajectory-sized storage without copying the complete accepted
        // history for each discarded heat root. Native stepping appends one row.
        let mut next = Self {
            model: self.model,
            state: SecondOrderState::new(self.state.t, &self.state.q, &self.state.v, &self.state.a),
            accepted: self.accepted.clone(),
            masses_kg: self.masses_kg.clone(),
        };
        let mut motion_evaluations = 0usize;
        let mut gas_trials = 0usize;
        let (thermal, step) = feedback::solve_scalar(
            exchange.heat_bracket_j,
            exchange.source_tolerance_j,
            exchange.maximum_iterations,
            |q, exchange_trials, cancelled| {
                let body_heat = body.internal_heat_j - q;
                let dh = body_heat / body.mass_kg;
                let h1 = h0 + dh;
                let gas_heat = heat_j + q;
                if ![body_heat, dh, h1, gas_heat].iter().all(|x| x.is_finite())
                    || (body_heat != 0.0 && (dh == 0.0 || h1 == h0))
                    || (q != 0.0 && (body_heat == body.internal_heat_j || gas_heat == heat_j))
                {
                    return Err(GasMembraneContactError::Invalid(
                        "heat exchange disappears or overflows in original energy coordinates",
                    ));
                }
                let thermal_midpoint = body
                    .curve
                    .state_at_specific_enthalpy(f64::midpoint(h0, h1))?;
                let endpoint = body.curve.state_at_specific_enthalpy(h1)?;
                let thermal_mean = body.curve.thermodynamic_mean_temperature(
                    h0,
                    h1,
                    body.curve.knots().len() - 1,
                    &mut *cancelled,
                )?;
                poll(cancelled)?;
                next.state.t = self.state.t;
                next.state.steps = self.state.steps;
                next.state.q.copy_from_slice(&self.state.q);
                next.state.v.copy_from_slice(&self.state.v);
                next.state.a.copy_from_slice(&self.state.a);
                next.state.history.clear();
                next.accepted = self.accepted.clone();
                let (motion, mean) = next.step_with_sample::<GasDiscreteGradientStep>(
                    duration_s,
                    dead_load_n,
                    gas_heat,
                    exterior_heat_j,
                    controls,
                    true,
                    &mut *cancelled,
                )?;
                let mean = mean.expect("requested original first-gas thermodynamic sample");
                motion_evaluations = motion_evaluations.checked_add(motion.evaluations).ok_or(
                    GasMembraneContactError::Invalid("contact evaluation count overflow"),
                )?;
                gas_trials = gas_trials.checked_add(motion.gas_trials).ok_or(
                    GasMembraneContactError::Invalid("contact gas trial count overflow"),
                )?;
                poll(cancelled)?;
                let sample =
                    resolve(&mean, thermal_midpoint).map_err(GasMembraneContactError::Interface)?;
                let g = sample.conductance_w_k;
                let tg = mean.temperature_k();
                let tb = thermal_mean.temperature_k;
                let delta = tb - tg;
                let power = g * delta;
                let evaluated = duration_s * power;
                let residual = q - evaluated;
                if !g.is_finite()
                    || g < 0.0
                    || ![power, evaluated, residual].iter().all(|x| x.is_finite())
                    || (g != 0.0 && delta != 0.0 && power == 0.0)
                    || (power != 0.0 && evaluated == 0.0)
                {
                    return Err(GasMembraneContactError::Invalid(
                        "contact conductance or heat is negative, nonfinite or unrepresentable",
                    ));
                }
                let stored = body.mass_kg * (endpoint.specific_enthalpy_j_kg() - h0);
                let body_residual = (Dd::from_f64(stored) + Dd::from_f64(q)
                    - Dd::from_f64(body.internal_heat_j))
                .to_f64();
                let other_storage = motion
                    .exterior_gas
                    .map_or(0.0, |gas| gas.stored_energy_change_j);
                let total = (Dd::from_f64(motion.kinetic_energy_change_j)
                    + Dd::from_f64(motion.membrane_energy_change_j)
                    + Dd::from_f64(motion.gas.stored_energy_change_j)
                    + Dd::from_f64(other_storage)
                    + Dd::from_f64(stored)
                    + Dd::from_f64(motion.exterior_pressure_work_j)
                    - Dd::from_f64(motion.dead_load_work_j)
                    - Dd::from_f64(heat_j)
                    - Dd::from_f64(exterior_heat_j.unwrap_or(0.0))
                    - Dd::from_f64(body.internal_heat_j))
                .to_f64();
                let body_entropy = body.mass_kg * thermal_mean.specific_entropy_change_j_kg_k;
                let other_entropy = motion
                    .exterior_gas
                    .map_or(0.0, |gas| gas.entropy_change_j_k);
                let matter_entropy = (Dd::from_f64(motion.gas.entropy_change_j_k)
                    + Dd::from_f64(other_entropy)
                    + Dd::from_f64(body_entropy))
                .to_f64();
                let gas_source_entropy = heat_j / tg;
                let body_source_entropy = body.internal_heat_j / tb;
                let (other_source_entropy, other_allowance) = match &self.accepted.exterior_gas {
                    Some(gas) => {
                        let temperature = motion
                            .exterior_sampling_temperature_k
                            .expect("sealed source sample");
                        (
                            exterior_heat_j.expect("admitted opposing heat") / temperature,
                            gas.mass_kg()
                                * controls
                                    .exterior_gas
                                    .expect("admitted opposing controls")
                                    .energy_tolerance_j_kg
                                / temperature,
                        )
                    }
                    None => (0.0, 0.0),
                };
                let prescribed_entropy = (Dd::from_f64(gas_source_entropy)
                    + Dd::from_f64(body_source_entropy)
                    + Dd::from_f64(other_source_entropy))
                .to_f64();
                // Stable nonnegative Fourier-contact production, evaluated from
                // the actual law rather than replacing it with the solved Q.
                let production = (evaluated / tg) * (delta / tb);
                let entropy_residual = (Dd::from_f64(matter_entropy)
                    - Dd::from_f64(prescribed_entropy)
                    - Dd::from_f64(production))
                .to_f64();
                let allowance = self.accepted.gas.mass_kg() * controls.gas.energy_tolerance_j_kg
                    / tg
                    + other_allowance
                    + exchange.energy_tolerance_j / tb
                    + exchange.source_tolerance_j * ((delta / tg) / tb).abs()
                    + 128.0
                        * f64::EPSILON
                        * (motion.gas.entropy_change_j_k.abs()
                            + other_entropy.abs()
                            + body_entropy.abs()
                            + gas_source_entropy.abs()
                            + body_source_entropy.abs()
                            + other_source_entropy.abs()
                            + production.abs());
                if ![
                    stored,
                    body_residual,
                    total,
                    body_entropy,
                    matter_entropy,
                    prescribed_entropy,
                    production,
                    entropy_residual,
                    allowance,
                ]
                .iter()
                .all(|x| x.is_finite())
                    || production < 0.0
                    || allowance < 0.0
                    || body_residual.abs() > exchange.energy_tolerance_j
                    || total.abs() > exchange.energy_tolerance_j
                {
                    return Err(GasMembraneContactError::Invalid(
                        "finite-body storage or combined energy/entropy account does not close",
                    ));
                }
                poll(cancelled)?;
                Ok((
                    (
                        UniformEnthalpyStep {
                            state: endpoint,
                            external_heat_j: -q,
                            energy_residual_tolerance_j: exchange.energy_tolerance_j,
                            report: sample.source,
                        },
                        GasMembraneContactStep {
                            motion,
                            gas_thermodynamic_mean: mean,
                            thermal_midpoint,
                            thermal_mean_temperature_k: tb,
                            conductance_w_k: g,
                            exchanged_heat_j: q,
                            evaluated_heat_j: evaluated,
                            source_residual_j: residual,
                            thermal_energy_change_j: stored,
                            thermal_energy_residual_j: body_residual,
                            total_energy_residual_j: total,
                            thermal_entropy_change_j_k: body_entropy,
                            total_entropy_change_j_k: matter_entropy,
                            prescribed_source_entropy_j_k: prescribed_entropy,
                            contact_entropy_production_j_k: production,
                            entropy_residual_j_k: entropy_residual,
                            entropy_tolerance_j_k: allowance,
                            exchange_trials,
                            motion_evaluations,
                            gas_trials,
                        },
                    ),
                    residual,
                ))
            },
            &mut cancelled,
        )
        .map_err(|error| match error {
            ScalarError::Evaluation(error) => error,
            ScalarError::Invalid(reason) => GasMembraneContactError::Invalid(reason),
            ScalarError::Cancelled => GasMembraneContactError::Cancelled,
            ScalarError::OutsideBracket { lower, upper } => {
                GasMembraneContactError::OutsideBracket {
                    lower_residual_j: lower,
                    upper_residual_j: upper,
                }
            }
            ScalarError::NotConverged {
                iterations,
                residual,
            } => GasMembraneContactError::NotConverged {
                iterations,
                residual_j: residual,
            },
        })?;
        if step.entropy_residual_j_k.abs() > step.entropy_tolerance_j_k
            || (heat_j == 0.0
                && body.internal_heat_j == 0.0
                && exterior_heat_j.unwrap_or(0.0) == 0.0
                && step.total_entropy_change_j_k < -step.entropy_tolerance_j_k)
        {
            return Err(GasMembraneContactError::Invalid(
                "contact entropy balance exceeds original energy-derived numerical allowance",
            ));
        }
        let entry = next
            .state
            .history
            .pop()
            .expect("one native accepted contact step");
        next.state.history.clone_from(&self.state.history);
        next.state.history.push(entry);
        poll(&mut cancelled)?;
        Ok(GasMembraneContactCandidate {
            motion: next,
            thermal,
            step,
        })
    }
}

#[cfg(test)]
mod tests;
