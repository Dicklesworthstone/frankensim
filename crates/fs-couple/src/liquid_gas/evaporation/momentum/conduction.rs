//! Mean-drift phase transfer, finite thermal storage and one native flow step.
//!
//! The shared conduction owner runs on the actual post-phase gas, including
//! receiver mixing heat. Its mass-preserving thermal stage retains that gas's
//! bulk momentum and kinetic inventory. No second phase, heat or flow solve is
//! introduced by attaching the mechanical account.

use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat,
    conduction::{
        FreeSurfaceBodyContact, FreeSurfaceBodyInterface, FreeSurfaceConductionControls,
        FreeSurfaceConductionFace, FreeSurfaceGasFace, FreeSurfaceRadiationFace,
        FreeSurfaceThermalBody,
    },
};
use fs_material::thermochemical_gas::energy::GasEnergySolve;

use super::super::{
    conduction::{ConductingLiquidGasReport, advance_native_with_momentum_and_conduction},
    multiple::DistributedLiquidGasSolve,
};
use super::{multiple::DistributedPhaseMomentumReport, *};

/// Original thermal/phase ledger plus its closed pre-flow momentum account.
#[derive(Clone, Debug)]
pub struct MovingGasConductingReport {
    /// Already includes the mechanical source defect and each receiver heat
    /// once. The full thermal flow defect must not be added a second time.
    pub native: ConductingLiquidGasReport,
    /// Ends after phase mixing, before the mass-preserving thermal stage and
    /// original pressure/flow step. This is not a whole-flow FSI account.
    pub momentum: DistributedPhaseMomentumReport,
}

/// Liquid, retained moving gas and all finite walls publish atomically.
pub struct MovingGasConductingCandidate<'s, 'm, 'catalog, 'b> {
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: MovingGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    bodies: Vec<FreeSurfaceThermalBody<'b>>,
    report: MovingGasConductingReport,
}

impl<'s, 'm, 'catalog, 'b> MovingGasConductingCandidate<'s, 'm, 'catalog, 'b> {
    pub const fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        &self.liquid
    }
    pub const fn gas(&self) -> &MovingGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub fn bodies(&self) -> &[FreeSurfaceThermalBody<'b>] {
        &self.bodies
    }
    pub const fn report(&self) -> &MovingGasConductingReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        MovingGasParcel,
        &'catalog ThermochemicalGasModel,
        Vec<FreeSurfaceThermalBody<'b>>,
        MovingGasConductingReport,
    ) {
        (
            self.liquid,
            self.gas,
            self.catalog,
            self.bodies,
            self.report,
        )
    }
}

/// Joint phase/momentum exchange followed by the full original thermal owner.
///
/// Liquid conduction, finite wall storage, wall gas films/radiation, attached
/// electrical heating and optional viscous heating use their original source
/// laws and admission budgets. Phase receiver mixing heats the original
/// receiving owner before thermal exchange. Thermal exchange preserves exact
/// gas mass, composition, volume and bulk momentum; it changes only internal
/// caloric storage. Optional pressure uses the final heated gas EOS during the
/// ONE original native flow step.
///
/// This remains a first-order phase/thermal/flow split on one native interval.
/// The homogeneous mean drift does not infer wall drag, radiation momentum,
/// compression work or a resolved gas jet. Pressure and stationary wall forces
/// retain their original owners and are outside the source momentum receipt.
/// Any source, conservation, geometry, speed or cancellation refusal discards
/// liquid, gas and finite walls together.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_surfaces_with_moving_gas_and_conduction<'s, 'm, 'catalog, 'b>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: &MovingGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: &[LiquidGasEvaporationSurface<'_>],
    bodies: &[FreeSurfaceThermalBody<'b>],
    faces: &[FreeSurfaceConductionFace],
    heat: &[FreeSurfaceCellHeat],
    contacts: &[FreeSurfaceBodyContact],
    interfaces: &[FreeSurfaceBodyInterface<'_>],
    exposed: &[FreeSurfaceGasFace],
    phase_controls: DistributedLiquidGasSolve,
    thermal_controls: FreeSurfaceConductionControls,
    gas_controls: GasEnergySolve,
    maximum_density_excursion: Option<f64>,
    with_viscous_heating: bool,
    gate: &CancelGate,
) -> Result<MovingGasConductingCandidate<'s, 'm, 'catalog, 'b>, EvaporatingLiquidGasError> {
    advance_liquid_surfaces_with_moving_gas_and_radiative_conduction(
        liquid,
        gas,
        catalog,
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
        gas_controls,
        maximum_density_excursion,
        with_viscous_heating,
        gate,
    )
}

/// Include original resolved-cell liquid radiation in the same transaction.
///
/// Radiation is external supply exactly once; gas films exchange internal
/// heat reciprocally. The original resolved thermal-boundary-layer or justified
/// lumped-temperature requirement, solved-surface-history refusals and all
/// limits of [`advance_liquid_surfaces_with_moving_gas_and_conduction`] apply.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_surfaces_with_moving_gas_and_radiative_conduction<'s, 'm, 'catalog, 'b>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: &MovingGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
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
    gas_controls: GasEnergySolve,
    maximum_density_excursion: Option<f64>,
    with_viscous_heating: bool,
    gate: &CancelGate,
) -> Result<MovingGasConductingCandidate<'s, 'm, 'catalog, 'b>, EvaporatingLiquidGasError> {
    let (candidate, motion) = advance_native_with_momentum_and_conduction(
        liquid,
        gas.parcel(),
        catalog,
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
        gas_controls,
        maximum_density_excursion,
        with_viscous_heating,
        Some(gas.momentum_kg_m_s()),
        gate,
    )?;
    let momentum = motion.ok_or(EvaporatingLiquidGasError::Invalid(
        "conducting moving gas omitted its phase momentum account",
    ))?;
    let (liquid, parcel, catalog, bodies, native) = candidate.into_parts();
    // Shared thermal admission requires bit-identical post-phase mass. Keeping
    // the exact post-phase momentum therefore keeps its kinetic inventory too.
    let gas = MovingGasParcel::from_momentum(parcel, momentum.gas_momentum_after_kg_m_s)?;
    poll(gate)?;
    Ok(MovingGasConductingCandidate {
        liquid,
        gas,
        catalog,
        bodies,
        report: MovingGasConductingReport { native, momentum },
    })
}
