//! Finite material heat exchange on the original elastic gas boundary.
//!
//! One signed Q enters the gas and leaves the finite enthalpy body. Every
//! heat trial repeats the existing gas/pressure/modal solve from its accepted
//! history. The original caloric, EOS, midpoint pressure, structural storage,
//! damping and source transport owners remain authoritative for ideal and
//! nonideal gases. No separate piston or heat-capacity surrogate is introduced.

use super::{GasValveMotion, ModalGasError, ModalGasSolve, ModalGasStep};
use crate::cavity_phs::piston::thermal::layers::{LayerConductionError, MaterialConductionLayers};
use crate::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve, feedback};
use crate::cavity_phs::piston::{PistonGas, PistonGasState};
use fs_exec::Cx;
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::{
    EquilibriumPhaseState, PhaseStateError, UniformEnthalpyStep, UniformEnthalpyStepInput,
};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

/// Preserve the original motion, caloric and interface refusal.
#[derive(Debug)]
pub enum ModalGasThermalError<E> {
    Motion(ModalGasError),
    Thermal(PhaseStateError),
    Interface(E),
    Invalid(&'static str),
    OutsideBracket {
        lower_residual_j: f64,
        upper_residual_j: f64,
    },
    NotConverged {
        iterations: usize,
        residual_j: f64,
    },
    Cancelled,
}

impl<E> From<ModalGasError> for ModalGasThermalError<E> {
    fn from(error: ModalGasError) -> Self {
        match error {
            ModalGasError::Cancelled => Self::Cancelled,
            other => Self::Motion(other),
        }
    }
}
impl<E> From<PhaseStateError> for ModalGasThermalError<E> {
    fn from(error: PhaseStateError) -> Self {
        Self::Thermal(error)
    }
}
impl<E: core::fmt::Display> core::fmt::Display for ModalGasThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Motion(error) => write!(f, "modal gas heat: {error}"),
            Self::Thermal(error) => write!(f, "modal gas thermal material: {error}"),
            Self::Interface(error) => write!(f, "modal gas heat interface: {error}"),
            Self::Invalid(message) => write!(f, "modal gas heat: {message}"),
            Self::OutsideBracket {
                lower_residual_j,
                upper_residual_j,
            } => write!(
                f,
                "modal gas heat not bracketed: {lower_residual_j}, {upper_residual_j} J"
            ),
            Self::NotConverged {
                iterations,
                residual_j,
            } => write!(
                f,
                "modal gas heat did not converge in {iterations} iterations: {residual_j} J"
            ),
            Self::Cancelled => write!(f, "modal gas heat cancelled"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ModalGasThermalError<E> {}

/// Actual gas, elastic and caloric work, with the one inter-body heat transfer.
#[derive(Clone, Debug)]
pub struct ModalGasThermalStep<S: PistonGasState = ResolvedThermochemicalGasState> {
    pub coupled: ModalGasStep<S>,
    /// State at the midpoint BODY ENTHALPY, including any liquid mass fraction.
    pub thermal_midpoint: EquilibriumPhaseState,
    pub conductance_w_k: f64,
    /// Signed heat into the gas and out of the body [J].
    pub exchanged_heat_j: f64,
    /// Q - dt G(T_body,mid - T_gas,mid) [J].
    pub source_residual_j: f64,
    pub body_energy_change_j: f64,
    pub body_energy_residual_j: f64,
    /// dU_gas + dE_modal + D_modal + dH_body + W_exterior - Q_gas - Q_body.
    /// Internal Q cancels; the original modal damping remains an external loss.
    pub total_energy_residual_j: f64,
    pub exchange_trials: usize,
    /// Complete native gas/elastic trial calls, summed over ALL heat probes.
    pub coupled_trials: usize,
}

/// Publish the original motion/gas and thermal endpoint together.
/// No accepted owner or clock is mutated while constructing this proposal.
#[derive(Clone, Debug)]
pub struct ModalGasThermalCandidate<S, G: PistonGas = ClosedGasParcel> {
    motion: GasValveMotion<G>,
    thermal: UniformEnthalpyStep<S>,
    step: ModalGasThermalStep<G::State>,
}
impl<S, G: PistonGas> ModalGasThermalCandidate<S, G> {
    pub const fn motion(&self) -> &GasValveMotion<G> {
        &self.motion
    }
    pub const fn thermal(&self) -> &UniformEnthalpyStep<S> {
        &self.thermal
    }
    pub const fn step(&self) -> &ModalGasThermalStep<G::State> {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        GasValveMotion<G>,
        UniformEnthalpyStep<S>,
        ModalGasThermalStep<G::State>,
    ) {
        (self.motion, self.thermal, self.step)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll<E>(cx: &Cx<'_>) -> Result<(), ModalGasThermalError<E>> {
    cx.checkpoint().map_err(|_| ModalGasThermalError::Cancelled)
}
fn feedback_error<E>(
    error: feedback::ScalarError<ModalGasThermalError<E>>,
) -> ModalGasThermalError<E> {
    match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Invalid(message) => ModalGasThermalError::Invalid(message),
        feedback::ScalarError::Cancelled => ModalGasThermalError::Cancelled,
        feedback::ScalarError::OutsideBracket { lower, upper } => {
            ModalGasThermalError::OutsideBracket {
                lower_residual_j: lower,
                upper_residual_j: upper,
            }
        }
        feedback::ScalarError::NotConverged {
            iterations,
            residual,
        } => ModalGasThermalError::NotConverged {
            iterations,
            residual_j: residual,
        },
    }
}

impl<G: PistonGas> GasValveMotion<G> {
    /// Joint finite-body heat, gas compression and native elastic motion.
    ///
    /// The two prescribed heats exclude the inter-body transfer Q. The
    /// repeatable interface callback sees the actual gas ENERGY/VOLUME midpoint
    /// and actual body ENTHALPY midpoint. It must not debit a live reservoir.
    /// A nonideal gas retains its original admitted EOS branch, calorics,
    /// density, compressibility and transport instead of an ideal reference.
    ///
    /// The finite body uses its original uniform-temperature, fixed-geometry,
    /// negligible-pV-work enthalpy convention. It is a separate thermal
    /// participant: heating it does not rebuild the structural basis, change
    /// stiffness, add thermal expansion or transfer history across phases.
    /// Native structural damping leaves the enclosure; it is not heat twice.
    ///
    /// Every Q probe restarts the accepted gas and modal state. Only a complete
    /// candidate whose source, body, pressure-work and enclosing energy gates
    /// pass is returned. Caller publication must adopt both endpoints together.
    /// All source-domain, trial-budget, timestep and cancellation gates remain
    /// active. Native midpoint coupling needs temporal refinement; neither
    /// equilibrium monotonicity nor an entropy theorem is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_candidate<S, E>(
        &self,
        gas_heat_j: f64,
        exterior_pressure_pa: f64,
        body: UniformEnthalpyStepInput<'_>,
        controls: ModalGasSolve,
        exchange: HeatExchangeSolve,
        mut resolve: impl FnMut(&G::State, EquilibriumPhaseState) -> Result<ConductanceSample<S>, E>,
        cx: &Cx<'_>,
    ) -> Result<ModalGasThermalCandidate<S, G>, ModalGasThermalError<E>> {
        poll(cx)?;
        let invalid = ModalGasThermalError::Invalid;
        feedback::admit::<E>(exchange)
            .map_err(|_| invalid("invalid heat bracket or work controls"))?;
        if !positive(body.mass_kg)
            || !positive(body.volume_m3)
            || !positive(body.surface_area_m2)
            || !positive(body.duration_s)
            || body.duration_s != self.motion().model().sample_period_s()
            || !body.internal_heat_j.is_finite()
            || !gas_heat_j.is_finite()
        {
            return Err(invalid(
                "finite heat, positive body geometry/mass and one native sample clock required",
            ));
        }
        let h0 = body.initial.specific_enthalpy_j_kg();
        if body.curve.state_at_specific_enthalpy(h0)? != body.initial {
            return Err(invalid(
                "finite body initial state belongs to a different material chart",
            ));
        }
        let mut coupled_trials = 0usize;
        let candidate = feedback::solve_scalar(
            exchange.heat_bracket_j, exchange.source_tolerance_j, exchange.maximum_iterations,
            |q, exchange_trials, _| {
                let body_heat = body.internal_heat_j - q;
                let gas_heat = gas_heat_j + q;
                let dh = body_heat / body.mass_kg;
                let h1 = h0 + dh;
                let hm = h0 + 0.5 * dh;
                if !body_heat.is_finite() || !gas_heat.is_finite() || !dh.is_finite()
                    || !h1.is_finite() || !hm.is_finite()
                    || (body_heat != 0.0 && (dh == 0.0 || h1 == h0))
                    || (q != 0.0 && (body_heat == body.internal_heat_j || gas_heat == gas_heat_j))
                {
                    return Err(invalid("reciprocal heat disappears or overflows in physical energy coordinates"));
                }
                let thermal_midpoint = body.curve.state_at_specific_enthalpy(hm)?;
                let endpoint = body.curve.state_at_specific_enthalpy(h1)?;
                let mut motion = self.clone();
                let coupled = motion.step(gas_heat, exterior_pressure_pa, controls, cx)?;
                coupled_trials = coupled_trials.checked_add(coupled.trials)
                    .ok_or_else(|| invalid("coupled heat trial count overflow"))?;
                poll(cx)?;
                let sample = resolve(&coupled.midpoint, thermal_midpoint)
                    .map_err(ModalGasThermalError::Interface)?;
                let g = sample.conductance_w_k;
                let delta = thermal_midpoint.temperature_k() - coupled.midpoint.temperature_k();
                let power = g * delta;
                let heat = body.duration_s * power;
                let residual = q - heat;
                if !g.is_finite() || g < 0.0 || !power.is_finite() || !heat.is_finite()
                    || !residual.is_finite() || (g != 0.0 && delta != 0.0 && power == 0.0)
                    || (power != 0.0 && heat == 0.0)
                {
                    return Err(invalid("source conductance or midpoint heat flux is not nonnegative and representable"));
                }
                let stored = body.mass_kg * (endpoint.specific_enthalpy_j_kg() - h0);
                let body_residual = stored + q - body.internal_heat_j;
                let total = coupled.gas.stored_energy_change_j + coupled.modal_energy_change_j
                    + coupled.motion.modal.viscous_dissipation_j + stored
                    + coupled.exterior_pressure_work_j - gas_heat_j - body.internal_heat_j;
                if [stored, body_residual, total].iter().any(|x| !x.is_finite())
                    || body_residual.abs() > exchange.energy_tolerance_j
                    || total.abs() > exchange.energy_tolerance_j
                {
                    return Err(invalid("finite-body or enclosing gas/elastic/thermal energy does not close"));
                }
                Ok((ModalGasThermalCandidate {
                    motion,
                    thermal: UniformEnthalpyStep { state: endpoint, external_heat_j: -q,
                        energy_residual_tolerance_j: exchange.energy_tolerance_j,
                        report: sample.source },
                    step: ModalGasThermalStep { coupled, thermal_midpoint, conductance_w_k: g,
                        exchanged_heat_j: q, source_residual_j: residual,
                        body_energy_change_j: stored, body_energy_residual_j: body_residual,
                        total_energy_residual_j: total, exchange_trials, coupled_trials },
                }, residual))
            },
            &mut || cx.checkpoint().is_err(),
        ).map_err(feedback_error)?;
        poll(cx)?;
        Ok(candidate)
    }

    /// Resolve the original series layers at the same two energy midpoints.
    /// The finite body's conductivity and enthalpy retain the same source card;
    /// the actual ideal or nonideal gas supplies its own transport state.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_layers_candidate(
        &self,
        gas_heat_j: f64,
        exterior_pressure_pa: f64,
        body: UniformEnthalpyStepInput<'_>,
        controls: ModalGasSolve,
        exchange: HeatExchangeSolve,
        layers: MaterialConductionLayers<'_>,
        cx: &Cx<'_>,
    ) -> Result<
        ModalGasThermalCandidate<ResolvedMaterialStatePoint, G>,
        ModalGasThermalError<LayerConductionError>,
    >
    where
        G::State: FluidTransport,
    {
        self.thermal_candidate(
            gas_heat_j,
            exterior_pressure_pa,
            body,
            controls,
            exchange,
            |gas, thermal| layers.resolve(gas, thermal),
            cx,
        )
    }
}
