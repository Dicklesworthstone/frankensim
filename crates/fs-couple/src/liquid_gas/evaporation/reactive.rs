//! Chemistry and phase transfer share the original finite gas and native liquid.
//!
//! The existing reaction owner first changes composition at conserved total
//! internal energy. The existing phase owner then exchanges that actual gas
//! with native liquid, optionally including reciprocal sensible heat, and
//! advances the original VOF flow exactly once. Species formation energy is
//! retained by its original owner and is never applied again as a heat source.

use crate::cavity_phs::piston::reactive::{
    GasReaction, GasReactionError, GasReactionSolve, GasReactionStep,
};
use fs_material::evaporation::gas::chamber::heat::GasLiquidChamberHeatSolve;

use super::multiple::{
    DistributedLiquidGasReport, DistributedLiquidGasSolve, LiquidGasHeatSurface,
    advance_liquid_surfaces_with_gas, advance_liquid_surfaces_with_gas_and_heat,
};
use super::{
    CancelGate, ClosedGasParcel, EvaporatingLiquidGasError, EvaporatingLiquidGasReport,
    EvaporatingLiquidGasSolve, GasFilmHeatSolve, LiquidGasEvaporationSurface,
    MaterialFreeSurfaceEnthalpy3, ThermochemicalGasModel, advance_liquid_with_gas,
    advance_liquid_with_gas_and_heat,
};

/// Chemistry and complete-interval limits, in addition to the original phase
/// and native controls. Stage errors cannot cancel to pass these budgets.
#[derive(Clone, Copy, Debug)]
pub struct ReactiveLiquidPhaseSolve {
    pub reaction: GasReactionSolve,
    /// Both chemical intermediate and final gas departures from initial T [K].
    pub maximum_gas_temperature_increment_k: f64,
    pub mass_tolerance_kg: f64,
    pub energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum ReactiveLiquidPhaseError {
    Reaction(GasReactionError),
    Phase(EvaporatingLiquidGasError),
    Invalid(&'static str),
    MassBalance {
        absolute_defect_kg: f64,
        tolerance_kg: f64,
    },
    EnergyBalance {
        absolute_defect_j: f64,
        tolerance_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for ReactiveLiquidPhaseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting gas/native-liquid phase transfer: {self:?}")
    }
}
impl std::error::Error for ReactiveLiquidPhaseError {}
impl From<GasReactionError> for ReactiveLiquidPhaseError {
    fn from(error: GasReactionError) -> Self {
        if error == GasReactionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Reaction(error)
        }
    }
}
impl From<EvaporatingLiquidGasError> for ReactiveLiquidPhaseError {
    fn from(error: EvaporatingLiquidGasError) -> Self {
        match error {
            EvaporatingLiquidGasError::Cancelled => Self::Cancelled,
            error => Self::Phase(error),
        }
    }
}

/// Original owner reports and the complete reaction/heat/phase/flow account.
#[derive(Clone, Debug)]
pub struct ReactiveLiquidPhaseReport<R> {
    /// One retained native timestep, also used by chemistry and phase transfer.
    pub duration_s: f64,
    pub reaction: GasReactionStep,
    /// Single- or multiple-surface report, beginning AFTER the reaction stage.
    pub phase: R,
    /// Original chemical mass defect plus actual phase/native transport defect.
    pub mass_residual_kg: f64,
    /// Original absolute phase/native defects plus absolute chemical mass error.
    pub absolute_mass_defect_kg: f64,
    /// Chemical stored-energy defect plus source-relative phase/native defect.
    /// Native external heat is already subtracted by its original flow owner.
    pub energy_residual_j: f64,
    /// Sum of absolute chemical, phase, native caloric and transport defects.
    pub absolute_energy_defect_j: f64,
}

/// No accepted input changes until the caller publishes both original owners.
pub struct ReactiveLiquidPhaseCandidate<'s, 'm, 'catalog, R> {
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    report: ReactiveLiquidPhaseReport<R>,
}
impl<'s, 'm, 'catalog, R> ReactiveLiquidPhaseCandidate<'s, 'm, 'catalog, R> {
    pub const fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        &self.liquid
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub const fn report(&self) -> &ReactiveLiquidPhaseReport<R> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        ReactiveLiquidPhaseReport<R>,
    ) {
        (self.liquid, self.gas, self.catalog, self.report)
    }
}

// Only the existing original reports implement this accounting seam. In
// particular, chemical production of the volatile species is not mistaken for
// a phase mass defect: the phase report begins at the reacted inventory.
trait PhaseAccount {
    fn duration_s(&self) -> f64;
    fn mass_residual_kg(&self) -> f64;
    fn absolute_mass_defect_kg(&self) -> f64;
    fn energy_residual_j(&self) -> f64;
    fn absolute_energy_defect_j(&self) -> f64;
}
impl PhaseAccount for EvaporatingLiquidGasReport {
    fn duration_s(&self) -> f64 {
        self.duration_s
    }
    fn mass_residual_kg(&self) -> f64 {
        self.mass_residual_kg
    }
    fn absolute_mass_defect_kg(&self) -> f64 {
        self.mass_residual_kg
            .abs()
            .max(self.volatile_mass_residual_kg.abs())
    }
    fn energy_residual_j(&self) -> f64 {
        self.energy_residual_j
    }
    fn absolute_energy_defect_j(&self) -> f64 {
        // The combined exchange receipt already sums heat and phase errors.
        // Keep their original absolute defects separate at this stricter gate.
        self.phase.energy_residual_j.abs()
            + self
                .heat
                .as_ref()
                .map_or(0.0, |heat| heat.energy_residual_j.abs())
            + self.native_phase_caloric_defect_j.abs()
            + self.flow.balance_defect_j.abs()
    }
}
impl PhaseAccount for DistributedLiquidGasReport {
    fn duration_s(&self) -> f64 {
        self.duration_s
    }
    fn mass_residual_kg(&self) -> f64 {
        self.mass_residual_kg
    }
    fn absolute_mass_defect_kg(&self) -> f64 {
        self.absolute_mass_defect_kg
    }
    fn energy_residual_j(&self) -> f64 {
        self.energy_residual_j
    }
    fn absolute_energy_defect_j(&self) -> f64 {
        self.absolute_energy_defect_j
    }
}

impl GasReaction {
    /// React the original finite gas, exchange phase mass/enthalpy on one
    /// native interface and advance its original flow once over the SAME dt.
    ///
    /// Chemistry uses this reaction's complete original source catalog; phase
    /// transfer uses that identical catalog and the actual reacted composition.
    /// Subsequent calls retain the remaining species, gas temperature and native
    /// populations. The existing source, geometry, inventory, caloric, native
    /// mass-representability and attached-owner admission gates remain active.
    /// Every failure discards the complete private candidate.
    ///
    /// This is a first-order reaction/phase/flow split. Chemical and phase
    /// stages share a physical interval; they do not advance two native steps.
    /// Refine that complete interval for accuracy. Gas remains well mixed in
    /// its explicitly fixed volume; no sealed-vessel displacement, gas traction,
    /// recoil, flame, nucleation or phase-front remap is inferred. The reaction
    /// entropy gate is local to its chemical stage, not a whole-system theorem.
    #[allow(clippy::too_many_arguments)]
    pub fn evaporating_liquid_candidate<'s, 'm>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &ClosedGasParcel,
        surface: LiquidGasEvaporationSurface<'_>,
        phase: EvaporatingLiquidGasSolve,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<
        ReactiveLiquidPhaseCandidate<'s, 'm, '_, EvaporatingLiquidGasReport>,
        ReactiveLiquidPhaseError,
    > {
        self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            advance_liquid_with_gas(liquid, reacted, self.catalog(), surface, phase, gate)
                .map(|candidate| candidate.into_parts())
        })
    }

    /// Same original chemistry, reciprocal sensible film heat, phase transfer
    /// and one native step. The supplied h excludes latent and vapor-carried
    /// sensible heat. All limits of [`Self::evaporating_liquid_candidate`] apply.
    #[allow(clippy::too_many_arguments)]
    pub fn evaporating_liquid_heat_candidate<'s, 'm>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &ClosedGasParcel,
        surface: LiquidGasEvaporationSurface<'_>,
        heat_transfer_coefficient_w_m2_k: f64,
        phase: EvaporatingLiquidGasSolve,
        heat: GasFilmHeatSolve,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<
        ReactiveLiquidPhaseCandidate<'s, 'm, '_, EvaporatingLiquidGasReport>,
        ReactiveLiquidPhaseError,
    > {
        self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            advance_liquid_with_gas_and_heat(
                liquid,
                reacted,
                self.catalog(),
                surface,
                heat_transfer_coefficient_w_m2_k,
                phase,
                heat,
                gate,
            )
            .map(|candidate| candidate.into_parts())
        })
    }

    /// React once, then exchange all original native patches with ONE common
    /// gas inventory before one flow step. No per-patch gas or liquid timeline
    /// is introduced. The original joint donor and canonical-order gates apply.
    #[allow(clippy::too_many_arguments)]
    pub fn evaporating_liquid_surfaces_candidate<'s, 'm>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &ClosedGasParcel,
        surfaces: &[LiquidGasEvaporationSurface<'_>],
        phase: DistributedLiquidGasSolve,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<
        ReactiveLiquidPhaseCandidate<'s, 'm, '_, DistributedLiquidGasReport>,
        ReactiveLiquidPhaseError,
    > {
        self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            advance_liquid_surfaces_with_gas(liquid, reacted, self.catalog(), surfaces, phase, gate)
                .map(|candidate| candidate.into_parts())
        })
    }

    /// React once, jointly exchange reciprocal sensible heat and phase mass
    /// with the common gas, then advance all original native cells once.
    /// The first-order split and scope of the single-surface entry remain.
    #[allow(clippy::too_many_arguments)]
    pub fn evaporating_liquid_surfaces_heat_candidate<'s, 'm>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &ClosedGasParcel,
        surfaces: &[LiquidGasHeatSurface<'_>],
        phase: DistributedLiquidGasSolve,
        heat: GasLiquidChamberHeatSolve,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<
        ReactiveLiquidPhaseCandidate<'s, 'm, '_, DistributedLiquidGasReport>,
        ReactiveLiquidPhaseError,
    > {
        self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            advance_liquid_surfaces_with_gas_and_heat(
                liquid,
                reacted,
                self.catalog(),
                surfaces,
                phase,
                heat,
                gate,
            )
            .map(|candidate| candidate.into_parts())
        })
    }

    fn native_phase_candidate<'s, 'm, 'catalog, R: PhaseAccount>(
        &'catalog self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &ClosedGasParcel,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
        exchange: impl FnOnce(
            &ClosedGasParcel,
        ) -> Result<
            (
                MaterialFreeSurfaceEnthalpy3<'s, 'm>,
                ClosedGasParcel,
                &'catalog ThermochemicalGasModel,
                R,
            ),
            EvaporatingLiquidGasError,
        >,
    ) -> Result<ReactiveLiquidPhaseCandidate<'s, 'm, 'catalog, R>, ReactiveLiquidPhaseError> {
        type Error = ReactiveLiquidPhaseError;
        let poll = || {
            if gate.is_requested() {
                Err(Error::Cancelled)
            } else {
                Ok(())
            }
        };
        poll()?;
        if !controls.maximum_gas_temperature_increment_k.is_finite()
            || controls.maximum_gas_temperature_increment_k <= 0.0
            || [controls.mass_tolerance_kg, controls.energy_tolerance_j]
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(Error::Invalid(
                "finite whole-interval mass/energy and positive gas temperature limits required",
            ));
        }
        let duration_s = liquid.flow().scaling().timestep_s();
        let initial_temperature = gas.state().state().temperature;
        let admit_temperature = |candidate: &ClosedGasParcel| {
            let increment = candidate.state().state().temperature - initial_temperature;
            if !increment.is_finite()
                || increment.abs() > controls.maximum_gas_temperature_increment_k
            {
                Err(Error::Invalid(
                    "chemical intermediate or phase endpoint exceeds the gas temperature allowance",
                ))
            } else {
                Ok(())
            }
        };
        let reaction =
            self.parcel_candidate(gas, duration_s, controls.reaction, || gate.is_requested())?;
        admit_temperature(reaction.gas())?;
        let (reacted, reaction) = reaction.into_parts();
        poll()?;
        let (next_liquid, next_gas, catalog, phase) = exchange(&reacted)?;
        admit_temperature(&next_gas)?;
        if phase.duration_s().to_bits() != duration_s.to_bits()
            || liquid.accepted_steps().checked_add(1) != Some(next_liquid.accepted_steps())
            || next_gas.volume_m3().to_bits() != gas.volume_m3().to_bits()
        {
            return Err(Error::Invalid(
                "chemistry/phase/flow lost its original single interval or fixed gas volume",
            ));
        }
        // Sum source-relative defects rather than subtracting large arbitrary
        // formation references again. Phase/native volatile balances start
        // after chemistry; stoichiometric production remains in reaction.
        let mass_residual_kg = reaction.mass_residual_kg + phase.mass_residual_kg();
        let absolute_mass_defect_kg =
            reaction.mass_residual_kg.abs() + phase.absolute_mass_defect_kg();
        let energy_residual_j = reaction.stored_energy_change_j + phase.energy_residual_j();
        let absolute_energy_defect_j =
            reaction.stored_energy_change_j.abs() + phase.absolute_energy_defect_j();
        if [
            mass_residual_kg,
            absolute_mass_defect_kg,
            energy_residual_j,
            absolute_energy_defect_j,
        ]
        .iter()
        .any(|value| !value.is_finite())
        {
            return Err(Error::Invalid(
                "combined reaction/phase/native inventory account is unrepresentable",
            ));
        }
        if absolute_mass_defect_kg > controls.mass_tolerance_kg {
            return Err(Error::MassBalance {
                absolute_defect_kg: absolute_mass_defect_kg,
                tolerance_kg: controls.mass_tolerance_kg,
            });
        }
        if absolute_energy_defect_j > controls.energy_tolerance_j {
            return Err(Error::EnergyBalance {
                absolute_defect_j: absolute_energy_defect_j,
                tolerance_j: controls.energy_tolerance_j,
            });
        }
        poll()?;
        Ok(ReactiveLiquidPhaseCandidate {
            liquid: next_liquid,
            gas: next_gas,
            catalog,
            report: ReactiveLiquidPhaseReport {
                duration_s,
                reaction,
                phase,
                mass_residual_kg,
                absolute_mass_defect_kg,
                energy_residual_j,
                absolute_energy_defect_j,
            },
        })
    }
}
