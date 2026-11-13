//! Chemistry and phase transfer share the original finite gas and native liquid.
//!
//! The existing reaction owner first changes composition at conserved total
//! internal energy. The existing phase owner then exchanges that actual gas
//! with native liquid, optionally including reciprocal sensible heat, and
//! advances the original VOF flow exactly once. The complete conduction entry
//! also retains finite walls, attached electrical heating and gas heat transfer.
//! Species formation energy stays with its original owner and is never applied
//! again as a heat source.

use super::conduction::{
    ConductingLiquidGasReport, advance_liquid_surfaces_with_gas_and_radiative_conduction,
};
use crate::cavity_phs::piston::reactive::{
    GasReaction, GasReactionError, GasReactionSolve, GasReactionStep,
};
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat,
    conduction::{
        FreeSurfaceBodyContact, FreeSurfaceBodyInterface, FreeSurfaceConductionControls,
        FreeSurfaceConductionFace, FreeSurfaceGasFace, FreeSurfaceRadiationFace,
        FreeSurfaceThermalBody,
    },
};
use fs_material::evaporation::gas::chamber::heat::GasLiquidChamberHeatSolve;

use super::multiple::{
    DistributedLiquidGasReport, DistributedLiquidGasSolve, LiquidGasHeatSurface,
    advance_liquid_surfaces_with_gas, advance_liquid_surfaces_with_gas_and_heat,
    advance_liquid_surfaces_with_pressurized_gas,
    advance_liquid_surfaces_with_pressurized_gas_and_heat,
};
use super::{
    CancelGate, ClosedGasParcel, EvaporatingLiquidGasError, EvaporatingLiquidGasReport,
    EvaporatingLiquidGasSolve, GasFilmHeatSolve, LiquidGasEvaporationSurface,
    MaterialFreeSurfaceEnthalpy3, ThermochemicalGasModel, advance_liquid_with_gas,
    advance_liquid_with_gas_and_heat, advance_liquid_with_pressurized_gas,
    advance_liquid_with_pressurized_gas_and_heat,
};

/// Chemistry and complete-interval limits, in addition to the original phase
/// and native controls. Stage errors cannot cancel to pass these budgets.
#[derive(Clone, Copy, Debug)]
pub struct ReactiveLiquidPhaseSolve {
    pub reaction: GasReactionSolve,
    /// Chemical/post-phase intermediates and final gas departures from initial T [K].
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
    /// Original phase/flow or phase/conduction/flow report, AFTER chemistry.
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
pub struct ReactiveLiquidPhaseCandidate<'s, 'm, 'catalog, R, G = ClosedGasParcel> {
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: G,
    catalog: &'catalog ThermochemicalGasModel,
    report: ReactiveLiquidPhaseReport<R>,
}
impl<'s, 'm, 'catalog, R, G> ReactiveLiquidPhaseCandidate<'s, 'm, 'catalog, R, G> {
    pub const fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        &self.liquid
    }
    pub const fn gas(&self) -> &G {
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
        G,
        &'catalog ThermochemicalGasModel,
        ReactiveLiquidPhaseReport<R>,
    ) {
        (self.liquid, self.gas, self.catalog, self.report)
    }
}

/// Chemical, phase and full native thermal/flow candidate with every finite
/// wall retained. All original owners publish together through this value.
pub struct ReactiveConductingLiquidGasCandidate<
    's,
    'm,
    'catalog,
    'b,
    R = ConductingLiquidGasReport,
    G = ClosedGasParcel,
> {
    pub(super) phase: ReactiveLiquidPhaseCandidate<'s, 'm, 'catalog, R, G>,
    pub(super) bodies: Vec<FreeSurfaceThermalBody<'b>>,
}
impl<'s, 'm, 'catalog, 'b, R, G> ReactiveConductingLiquidGasCandidate<'s, 'm, 'catalog, 'b, R, G> {
    pub fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        self.phase.liquid()
    }
    pub fn gas(&self) -> &G {
        self.phase.gas()
    }
    pub fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.phase.catalog()
    }
    pub fn bodies(&self) -> &[FreeSurfaceThermalBody<'b>] {
        &self.bodies
    }
    pub fn report(&self) -> &ReactiveLiquidPhaseReport<R> {
        self.phase.report()
    }
    pub fn into_parts(
        self,
    ) -> (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        G,
        &'catalog ThermochemicalGasModel,
        Vec<FreeSurfaceThermalBody<'b>>,
        ReactiveLiquidPhaseReport<R>,
    ) {
        let (liquid, gas, catalog, report) = self.phase.into_parts();
        (liquid, gas, catalog, self.bodies, report)
    }
}

// Chemistry still runs exclusively on the original parcel owner. A gas with
// an explicit mechanical inventory must reattach it only after that owner's
// mass-preserving reaction candidate has passed every original admission.
pub(super) trait PhaseGas: Sized {
    fn parcel(&self) -> &ClosedGasParcel;
    fn with_reacted(&self, reacted: ClosedGasParcel) -> Result<Self, ReactiveLiquidPhaseError>;
}
impl PhaseGas for ClosedGasParcel {
    fn parcel(&self) -> &ClosedGasParcel {
        self
    }
    fn with_reacted(&self, reacted: ClosedGasParcel) -> Result<Self, ReactiveLiquidPhaseError> {
        Ok(reacted)
    }
}

// Only the existing original reports implement this accounting seam. In
// particular, chemical production of the volatile species is not mistaken for
// a phase mass defect: the phase report begins at the reacted inventory.
pub(super) trait PhaseAccount {
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

impl PhaseAccount for ConductingLiquidGasReport {
    fn duration_s(&self) -> f64 {
        self.phase.duration_s
    }
    fn mass_residual_kg(&self) -> f64 {
        self.phase.mass_residual_kg
    }
    fn absolute_mass_defect_kg(&self) -> f64 {
        self.phase.absolute_mass_defect_kg
    }
    fn energy_residual_j(&self) -> f64 {
        self.energy_residual_j
    }
    fn absolute_energy_defect_j(&self) -> f64 {
        self.absolute_energy_defect_j
    }
}

impl<R> PhaseAccount for ReactiveLiquidPhaseReport<R> {
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

    /// React, exchange phase mass and energy, then apply the resulting gas
    /// pressure during the SAME native flow step.
    ///
    /// The original chemical and phase owners determine the final finite-gas
    /// composition, mass and temperature. Its actual EOS pressure drives the
    /// existing native missing-link reconstruction; neither initial nor
    /// chemistry-only pressure is substituted. The phase pressure receipt
    /// begins after chemistry and records the pressure actually applied.
    ///
    /// Every exposed native gas link must share this one well-mixed, explicitly
    /// fixed-volume parcel. The original pressure owner requires a fixed
    /// pressure reference, stationary geometry and a density excursion in
    /// (0,1). Phase transfer retains the original mass/enthalpy account and
    /// single native clock. Pressure loading adds no inferred caloric heat;
    /// this fixed-volume split does not model sealed-vessel displacement,
    /// gas compression work or a combined mechanical-energy theorem.
    #[allow(clippy::too_many_arguments)]
    pub fn pressurized_evaporating_liquid_candidate<'s, 'm>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &ClosedGasParcel,
        surface: LiquidGasEvaporationSurface<'_>,
        maximum_density_excursion: f64,
        phase: EvaporatingLiquidGasSolve,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<
        ReactiveLiquidPhaseCandidate<'s, 'm, '_, EvaporatingLiquidGasReport>,
        ReactiveLiquidPhaseError,
    > {
        self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            advance_liquid_with_pressurized_gas(
                liquid,
                reacted,
                self.catalog(),
                surface,
                maximum_density_excursion,
                phase,
                gate,
            )
            .map(|candidate| candidate.into_parts())
        })
    }

    /// Include reciprocal sensible heat in the original phase transaction,
    /// then load native flow with the post-heat/post-phase gas EOS pressure.
    /// All scope and admission limits of
    /// [`Self::pressurized_evaporating_liquid_candidate`] remain in force.
    #[allow(clippy::too_many_arguments)]
    pub fn pressurized_evaporating_liquid_heat_candidate<'s, 'm>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &ClosedGasParcel,
        surface: LiquidGasEvaporationSurface<'_>,
        heat_transfer_coefficient_w_m2_k: f64,
        maximum_density_excursion: f64,
        phase: EvaporatingLiquidGasSolve,
        heat: GasFilmHeatSolve,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<
        ReactiveLiquidPhaseCandidate<'s, 'm, '_, EvaporatingLiquidGasReport>,
        ReactiveLiquidPhaseError,
    > {
        self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            advance_liquid_with_pressurized_gas_and_heat(
                liquid,
                reacted,
                self.catalog(),
                surface,
                heat_transfer_coefficient_w_m2_k,
                maximum_density_excursion,
                phase,
                heat,
                gate,
            )
            .map(|candidate| candidate.into_parts())
        })
    }

    /// React one common gas, exchange all original native phase surfaces,
    /// and use the one joint endpoint pressure for one native flow step.
    /// Joint donor limits and canonical surface ordering remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn pressurized_evaporating_liquid_surfaces_candidate<'s, 'm>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &ClosedGasParcel,
        surfaces: &[LiquidGasEvaporationSurface<'_>],
        maximum_density_excursion: f64,
        phase: DistributedLiquidGasSolve,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<
        ReactiveLiquidPhaseCandidate<'s, 'm, '_, DistributedLiquidGasReport>,
        ReactiveLiquidPhaseError,
    > {
        self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            advance_liquid_surfaces_with_pressurized_gas(
                liquid,
                reacted,
                self.catalog(),
                surfaces,
                maximum_density_excursion,
                phase,
                gate,
            )
            .map(|candidate| candidate.into_parts())
        })
    }

    /// Compose chemistry, joint reciprocal heat and phase transfer, then load
    /// every exposed native gas link with the SAME final gas pressure.
    /// Each original owner shares the retained interval and publishes only
    /// through the complete candidate; per-patch clocks are never advanced.
    #[allow(clippy::too_many_arguments)]
    pub fn pressurized_evaporating_liquid_surfaces_heat_candidate<'s, 'm>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &ClosedGasParcel,
        surfaces: &[LiquidGasHeatSurface<'_>],
        maximum_density_excursion: f64,
        phase: DistributedLiquidGasSolve,
        heat: GasLiquidChamberHeatSolve,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<
        ReactiveLiquidPhaseCandidate<'s, 'm, '_, DistributedLiquidGasReport>,
        ReactiveLiquidPhaseError,
    > {
        self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            advance_liquid_surfaces_with_pressurized_gas_and_heat(
                liquid,
                reacted,
                self.catalog(),
                surfaces,
                maximum_density_excursion,
                phase,
                heat,
                gate,
            )
            .map(|candidate| candidate.into_parts())
        })
    }

    /// React the original gas, exchange phase mass with all native patches,
    /// then run the original gas/liquid/wall conduction and native flow once.
    ///
    /// The full thermal owner retains wall radiation, attached electrical
    /// heating, optional molecular heat and optional post-thermal gas pressure.
    /// Those sources are not synthesized or applied by this composition layer.
    /// Chemistry, phase and thermal exchange use the SAME retained interval;
    /// the existing shared reaction account checks the complete candidate.
    ///
    /// Chemical and post-phase intermediate temperatures, plus the final gas
    /// temperature, obey the declared departure bound. Source, pressure,
    /// geometry, donor, caloric and attached-owner admission remain with their
    /// original implementations. Failures preserve liquid, gas and all walls.
    /// This is a first-order chemistry/phase/conduction/flow split with fixed
    /// gas volume; it does not infer gas compression work or volume displacement.
    #[allow(clippy::too_many_arguments)]
    pub fn conducting_evaporating_liquid_candidate<'s, 'm, 'b>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &ClosedGasParcel,
        surfaces: &[LiquidGasEvaporationSurface<'_>],
        bodies: &[FreeSurfaceThermalBody<'b>],
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        exposed: &[FreeSurfaceGasFace],
        phase_controls: DistributedLiquidGasSolve,
        thermal_controls: FreeSurfaceConductionControls,
        maximum_density_excursion: Option<f64>,
        with_viscous_heating: bool,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<ReactiveConductingLiquidGasCandidate<'s, 'm, '_, 'b>, ReactiveLiquidPhaseError>
    {
        self.radiative_conducting_evaporating_liquid_candidate(
            liquid,
            gas,
            surfaces,
            bodies,
            faces,
            heat,
            contacts,
            interfaces,
            exposed,
            &[],
            phase_controls,
            thermal_controls,
            maximum_density_excursion,
            with_viscous_heating,
            controls,
            gate,
        )
    }

    /// Include original resolved-cell-temperature liquid radiation in the
    /// complete chemistry/phase/conduction transaction.
    ///
    /// This uses the existing gray-surface radiation law on the native thermal
    /// owner. It neither supplies a second sensible-film coefficient nor changes
    /// the phase owner's bulk-caloric temperature convention. Solved-interface
    /// latent feedback and attached surface/contact-history refusals remain
    /// unchanged. Radiation is external supply once; all finite participants,
    /// original sources, budgets and the single native clock remain shared.
    #[allow(clippy::too_many_arguments)]
    pub fn radiative_conducting_evaporating_liquid_candidate<'s, 'm, 'b>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &ClosedGasParcel,
        surfaces: &[LiquidGasEvaporationSurface<'_>],
        bodies: &[FreeSurfaceThermalBody<'b>],
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        exposed: &[FreeSurfaceGasFace],
        radiation: &[FreeSurfaceRadiationFace],
        phase_controls: DistributedLiquidGasSolve,
        thermal_controls: FreeSurfaceConductionControls,
        maximum_density_excursion: Option<f64>,
        with_viscous_heating: bool,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<ReactiveConductingLiquidGasCandidate<'s, 'm, '_, 'b>, ReactiveLiquidPhaseError>
    {
        let mut accepted_bodies = None;
        let phase = self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            let candidate = advance_liquid_surfaces_with_gas_and_radiative_conduction(
                liquid,
                reacted,
                self.catalog(),
                surfaces,
                bodies,
                faces,
                heat,
                contacts,
                interfaces,
                exposed,
                radiation,
                phase_controls,
                thermal_controls,
                controls.reaction.gas,
                maximum_density_excursion,
                with_viscous_heating,
                gate,
            )?;
            let (liquid, next_gas, catalog, next_bodies, report) = candidate.into_parts();
            let intermediate =
                report.thermal_gas_temperature_k[0] - gas.state().state().temperature;
            if !intermediate.is_finite()
                || intermediate.abs() > controls.maximum_gas_temperature_increment_k
            {
                return Err(EvaporatingLiquidGasError::Invalid(
                    "post-phase gas exceeds the reactive temperature-departure allowance",
                ));
            }
            accepted_bodies = Some(next_bodies);
            Ok((liquid, next_gas, catalog, report))
        })?;
        let bodies = accepted_bodies.ok_or(ReactiveLiquidPhaseError::Invalid(
            "complete reactive phase/conduction candidate lost its finite walls",
        ))?;
        Ok(ReactiveConductingLiquidGasCandidate { phase, bodies })
    }

    pub(super) fn native_phase_candidate<'s, 'm, 'catalog, R: PhaseAccount, G: PhaseGas>(
        &'catalog self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &G,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
        exchange: impl FnOnce(
            &G,
        ) -> Result<
            (
                MaterialFreeSurfaceEnthalpy3<'s, 'm>,
                G,
                &'catalog ThermochemicalGasModel,
                R,
            ),
            EvaporatingLiquidGasError,
        >,
    ) -> Result<ReactiveLiquidPhaseCandidate<'s, 'm, 'catalog, R, G>, ReactiveLiquidPhaseError>
    {
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
        let initial_temperature = gas.parcel().state().state().temperature;
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
            self.parcel_candidate(gas.parcel(), duration_s, controls.reaction, || {
                gate.is_requested()
            })?;
        admit_temperature(reaction.gas())?;
        let (reacted, reaction) = reaction.into_parts();
        let reacted = gas.with_reacted(reacted)?;
        poll()?;
        let (next_liquid, next_gas, catalog, phase) = exchange(&reacted)?;
        admit_temperature(next_gas.parcel())?;
        if phase.duration_s().to_bits() != duration_s.to_bits()
            || liquid.accepted_steps().checked_add(1) != Some(next_liquid.accepted_steps())
            || next_gas.parcel().volume_m3().to_bits() != gas.parcel().volume_m3().to_bits()
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
