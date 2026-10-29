//! Endpoint gas chemistry feeds the original continuing material-string sound.
//!
//! The native string step supplies actual aerodynamic heat. The existing
//! implicit reaction owner resolves composition and temperature with that heat
//! on the same sample interval, then the original finite-gas publisher rebuilds
//! drag and compact pressure from the accepted mixture. No extra chemical heat
//! is added: species formation energy is already part of gas internal energy.

use super::{
    FiniteGasStringCandidate, FiniteGasStringError, FiniteGasStringFrame, FiniteGasStringHeat,
    FiniteGasStringSolve, ThermalMaterialStringRuntime,
};
use crate::cavity_phs::piston::reactive::{
    GasReaction, GasReactionSolve, heat::HeatedGasReactionStep,
};
use fs_material::thermochemical_gas::energy::{IsochoricGasParcel, parcel::ClosedGasParcel};

/// Original chemical/caloric work limits and the mechanical split's bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveFiniteGasStringSolve {
    /// Uses the native implicit mechanism's source, amount, entropy and work gates.
    pub reaction: GasReactionSolve,
    /// Maximum absolute final gas temperature change in this sample [K].
    /// This bounds the frozen-medium split; it is not a temporal error estimate.
    pub maximum_gas_temperature_increment_k: f64,
}

/// One accepted sample with the complete original gas/chemistry accounts.
#[derive(Debug)]
pub struct ReactiveFiniteGasStringFrame {
    /// Measured mechanical, solid and formation-inclusive gas energy balance.
    pub coupled: FiniteGasStringFrame,
    /// Endpoint species, rates, nonlinear residuals, heat and entropy accounts.
    pub chemistry: HeatedGasReactionStep,
}

impl ThermalMaterialStringRuntime {
    /// Advance one material-string sample in a reacting finite atmosphere.
    ///
    /// The mechanical sample freezes the accepted initial solid and gas laws,
    /// exactly as in `step_with_isochoric_gas`. Native aerodynamic dissipation
    /// and prescribed gas heat enter ONE implicit chemical/caloric endpoint at
    /// that same sample duration. The mechanism must be configured with
    /// `GasReaction::with_implicit_mechanism`; all original rate scales, source
    /// selection, species limits, residual tolerances and cancellation remain
    /// authoritative. There is no reaction-only substep or inferred mechanism.
    ///
    /// The original finite-gas publisher independently measures complete energy,
    /// rebuilds drag and diagnostic pressure using the final mixture, retains
    /// prior acceleration and advances one sample/epoch. Both accepted owners
    /// remain unchanged if chemistry, source, energy, observer or work admission
    /// fails. Retrying therefore cannot deposit heat or emit a sample twice.
    ///
    /// This is a first-order mechanical/thermal split with simultaneous heated
    /// chemistry in a well-mixed fixed-volume gas. It omits gas motion, flame or
    /// chamber propagation, pressure work, acoustic radiation back-reaction and
    /// reciprocal conductive exchange. Its compact microphone is diagnostic;
    /// neither its pressure nor chemical formation release is another heat input.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_reacting_isochoric_gas(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        reaction: &GasReaction,
        controls: ReactiveFiniteGasStringSolve,
    ) -> Result<ReactiveFiniteGasStringFrame, FiniteGasStringError> {
        let finite = FiniteGasStringSolve {
            gas: controls.reaction.gas,
            maximum_gas_temperature_increment_k: controls.maximum_gas_temperature_increment_k,
        };
        self.admit_finite_gas_step(cx, expected_epoch, gas, heat, finite)?;
        let mut candidate = self.clone();
        let string = candidate.step(cx, expected_epoch, forces, heat.string_heat_j)?;
        let drag_heat_j = string.vibration.dissipation.air_loss_j;
        let total_gas_heat = heat.gas_heat_j + drag_heat_j;
        if !total_gas_heat.is_finite()
            || (drag_heat_j != 0.0 && total_gas_heat == heat.gas_heat_j)
            || (heat.gas_heat_j != 0.0 && total_gas_heat == drag_heat_j)
        {
            return Err(FiniteGasStringError::Admission(
                "gas source and aerodynamic heat overflow or disappear when combined",
            ));
        }
        let (next_gas, chemistry) = reaction
            .heated_implicit_candidate(
                &ClosedGasParcel::from_isochoric(gas.clone()),
                self.mechanical.member.model.sample_period_s(),
                total_gas_heat,
                controls.reaction,
                || cx.checkpoint().is_err(),
            )?
            .into_parts();
        let coupled = self.publish_finite_gas(
            cx,
            forces,
            gas,
            heat,
            finite,
            FiniteGasStringCandidate {
                runtime: candidate,
                gas: next_gas.into_isochoric(),
                string,
                gas_step: chemistry.heat.clone(),
            },
        )?;
        Ok(ReactiveFiniteGasStringFrame { coupled, chemistry })
    }
}
