//! Original chemistry, retained gas momentum and native phase/thermal flow.
//!
//! Homogeneous chemistry changes the original composition and temperature at
//! conserved internal energy and exact parcel mass. It exerts no external
//! impulse, so bulk momentum and kinetic energy stay unchanged. The reacted
//! parcel then enters the same phase donor/mixing owner as a nonreacting gas;
//! chemical formation energy is never reapplied as heat or drift energy.

use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat,
    conduction::{
        FreeSurfaceBodyContact, FreeSurfaceBodyInterface, FreeSurfaceConductionControls,
        FreeSurfaceConductionFace, FreeSurfaceGasFace, FreeSurfaceRadiationFace,
        FreeSurfaceThermalBody,
    },
};
use fs_material::evaporation::gas::chamber::heat::GasLiquidChamberHeatSolve;

use super::super::{
    multiple::{DistributedLiquidGasSolve, LiquidGasHeatSurface},
    reactive::{
        PhaseAccount, PhaseGas, ReactiveConductingLiquidGasCandidate, ReactiveLiquidPhaseCandidate,
        ReactiveLiquidPhaseError, ReactiveLiquidPhaseSolve,
    },
};
use super::{
    conduction::{
        MovingGasConductingReport, advance_liquid_surfaces_with_moving_gas_and_radiative_conduction,
    },
    multiple::{
        MovingGasDistributedEvaporationReport, advance_liquid_surfaces_with_moving_gas,
        advance_liquid_surfaces_with_moving_gas_and_heat,
    },
    *,
};
use crate::cavity_phs::piston::reactive::GasReaction;

/// Shared atomic reaction/phase candidate retaining the original moving gas.
pub type ReactiveMovingGasCandidate<'s, 'm, 'catalog, R> =
    ReactiveLiquidPhaseCandidate<'s, 'm, 'catalog, R, MovingGasParcel>;

/// Shared atomic reaction/phase/thermal candidate retaining every finite wall.
pub type ReactiveConductingMovingGasCandidate<'s, 'm, 'catalog, 'b> =
    ReactiveConductingLiquidGasCandidate<
        's,
        'm,
        'catalog,
        'b,
        MovingGasConductingReport,
        MovingGasParcel,
    >;

impl PhaseGas for MovingGasParcel {
    fn parcel(&self) -> &ClosedGasParcel {
        self.parcel()
    }
    fn with_reacted(&self, reacted: ClosedGasParcel) -> Result<Self, ReactiveLiquidPhaseError> {
        // GasReaction already enforces this invariant on each accepted
        // chemical endpoint, including implicit/subcycled mechanisms. Require
        // it explicitly here before preserving an independent kinetic account.
        if reacted.mass_kg().to_bits() != self.parcel().mass_kg().to_bits()
            || reacted.volume_m3().to_bits() != self.parcel().volume_m3().to_bits()
        {
            return Err(ReactiveLiquidPhaseError::Invalid(
                "reaction changed moving gas mass or its fixed volume",
            ));
        }
        Ok(Self::from_momentum(reacted, self.momentum_kg_m_s())?)
    }
}

impl PhaseAccount for MovingGasEvaporationReport {
    fn duration_s(&self) -> f64 {
        self.native.duration_s()
    }
    fn mass_residual_kg(&self) -> f64 {
        self.native.mass_residual_kg()
    }
    fn absolute_mass_defect_kg(&self) -> f64 {
        self.native.absolute_mass_defect_kg()
    }
    fn energy_residual_j(&self) -> f64 {
        self.native.energy_residual_j()
    }
    fn absolute_energy_defect_j(&self) -> f64 {
        self.native.absolute_energy_defect_j() + self.momentum.energy_residual_j.abs()
    }
}
impl PhaseAccount for MovingGasDistributedEvaporationReport {
    fn duration_s(&self) -> f64 {
        self.native.duration_s()
    }
    fn mass_residual_kg(&self) -> f64 {
        self.native.mass_residual_kg()
    }
    fn absolute_mass_defect_kg(&self) -> f64 {
        self.native.absolute_mass_defect_kg()
    }
    fn energy_residual_j(&self) -> f64 {
        self.native.energy_residual_j()
    }
    fn absolute_energy_defect_j(&self) -> f64 {
        // Distributed native accounting already includes the mechanical defect.
        self.native.absolute_energy_defect_j()
    }
}
impl PhaseAccount for MovingGasConductingReport {
    fn duration_s(&self) -> f64 {
        self.native.duration_s()
    }
    fn mass_residual_kg(&self) -> f64 {
        self.native.mass_residual_kg()
    }
    fn absolute_mass_defect_kg(&self) -> f64 {
        self.native.absolute_mass_defect_kg()
    }
    fn energy_residual_j(&self) -> f64 {
        self.native.energy_residual_j()
    }
    fn absolute_energy_defect_j(&self) -> f64 {
        self.native.absolute_energy_defect_j()
    }
}

impl GasReaction {
    /// React the actual finite gas, then exchange phase and reciprocal drift.
    ///
    /// Optional sensible heat and optional native pressure use the original
    /// moving-phase entry. Chemistry, phase and flow share ONE native interval;
    /// the temperature-departure bound includes chemical and final mixed gas
    /// endpoints relative to the original gas. Original reaction integrators,
    /// donor/source/geometry admission and absolute defect budgets are retained.
    /// Any failure leaves both accepted owners and their clock untouched.
    #[allow(clippy::too_many_arguments)]
    pub fn moving_gas_evaporating_liquid_candidate<'s, 'm>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &MovingGasParcel,
        surface: LiquidGasEvaporationSurface<'_>,
        phase_controls: EvaporatingLiquidGasSolve,
        heat_controls: Option<(f64, GasFilmHeatSolve)>,
        maximum_density_excursion: Option<f64>,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<
        ReactiveMovingGasCandidate<'s, 'm, '_, MovingGasEvaporationReport>,
        ReactiveLiquidPhaseError,
    > {
        self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            advance_liquid_with_moving_gas(
                liquid,
                reacted,
                self.catalog(),
                surface,
                phase_controls,
                heat_controls,
                maximum_density_excursion,
                gate,
            )
            .map(MovingGasEvaporationCandidate::into_parts)
        })
    }

    /// One reacted donor inventory shared by every native phase interface.
    ///
    /// Gross condensation uses the original reacted gas drift; each evaporating
    /// donor uses its actual native liquid velocity. Retained and incoming gas
    /// mixes once before optional final-EOS pressure and the single flow step.
    #[allow(clippy::too_many_arguments)]
    pub fn moving_gas_evaporating_liquid_surfaces_candidate<'s, 'm>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &MovingGasParcel,
        surfaces: &[LiquidGasEvaporationSurface<'_>],
        phase_controls: DistributedLiquidGasSolve,
        maximum_density_excursion: Option<f64>,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<
        ReactiveMovingGasCandidate<'s, 'm, '_, MovingGasDistributedEvaporationReport>,
        ReactiveLiquidPhaseError,
    > {
        self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            advance_liquid_surfaces_with_moving_gas(
                liquid,
                reacted,
                self.catalog(),
                surfaces,
                phase_controls,
                maximum_density_excursion,
                gate,
            )
            .map(multiple::MovingGasDistributedEvaporationCandidate::into_parts)
        })
    }

    /// Include the original reciprocal sensible-film solve before joint phase
    /// transport. Its internal heat, receiver mixing and chemical formation
    /// energy retain separate original accounts and are each applied once.
    #[allow(clippy::too_many_arguments)]
    pub fn moving_gas_evaporating_liquid_surfaces_heat_candidate<'s, 'm>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &MovingGasParcel,
        surfaces: &[LiquidGasHeatSurface<'_>],
        phase_controls: DistributedLiquidGasSolve,
        heat_controls: GasLiquidChamberHeatSolve,
        maximum_density_excursion: Option<f64>,
        controls: ReactiveLiquidPhaseSolve,
        gate: &CancelGate,
    ) -> Result<
        ReactiveMovingGasCandidate<'s, 'm, '_, MovingGasDistributedEvaporationReport>,
        ReactiveLiquidPhaseError,
    > {
        self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            advance_liquid_surfaces_with_moving_gas_and_heat(
                liquid,
                reacted,
                self.catalog(),
                surfaces,
                phase_controls,
                heat_controls,
                maximum_density_excursion,
                gate,
            )
            .map(multiple::MovingGasDistributedEvaporationCandidate::into_parts)
        })
    }

    /// Chemistry, phase drift and the original finite gas/wall thermal timeline.
    /// All finite bodies publish with the liquid and retained moving gas.
    #[allow(clippy::too_many_arguments)]
    pub fn moving_gas_conducting_evaporating_liquid_candidate<'s, 'm, 'b>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &MovingGasParcel,
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
    ) -> Result<ReactiveConductingMovingGasCandidate<'s, 'm, '_, 'b>, ReactiveLiquidPhaseError>
    {
        self.moving_gas_radiative_conducting_evaporating_liquid_candidate(
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

    /// Add original resolved-cell liquid radiation to the complete transaction.
    ///
    /// The exact reacted gas enters joint phase/mixing, then the original
    /// finite-wall/gas conduction and optional pressure/flow owner. Chemical,
    /// post-phase/mixing and final thermal gas endpoints obey the departure
    /// allowance from the original temperature. No second flow step, formation
    /// heat or kinetic-to-internal correction is inserted into chemistry or the
    /// mass-preserving thermal stage.
    ///
    /// Original fixed-volume, gray-radiation, source, wall geometry and attached
    /// history limits remain. This first-order split retains a source-only
    /// phase momentum account, not a closed mechanical balance for pressure,
    /// radiation, wall forces or the entire native flow.
    #[allow(clippy::too_many_arguments)]
    pub fn moving_gas_radiative_conducting_evaporating_liquid_candidate<'s, 'm, 'b>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &MovingGasParcel,
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
    ) -> Result<ReactiveConductingMovingGasCandidate<'s, 'm, '_, 'b>, ReactiveLiquidPhaseError>
    {
        let mut accepted_bodies = None;
        let phase = self.native_phase_candidate(liquid, gas, controls, gate, |reacted| {
            let candidate = advance_liquid_surfaces_with_moving_gas_and_radiative_conduction(
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
            let intermediate = report.native.thermal_gas_temperature_k[0]
                - gas.parcel().state().state().temperature;
            if !intermediate.is_finite()
                || intermediate.abs() > controls.maximum_gas_temperature_increment_k
            {
                return Err(EvaporatingLiquidGasError::Invalid(
                    "post-phase moving gas exceeds the reactive temperature-departure allowance",
                ));
            }
            accepted_bodies = Some(next_bodies);
            Ok((liquid, next_gas, catalog, report))
        })?;
        let bodies = accepted_bodies.ok_or(ReactiveLiquidPhaseError::Invalid(
            "reactive moving gas phase/conduction candidate lost its finite walls",
        ))?;
        Ok(ReactiveConductingLiquidGasCandidate { phase, bodies })
    }
}
