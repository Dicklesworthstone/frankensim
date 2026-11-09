//! Finite-gas transport feeding distinct local phase/heat atmospheres.
//!
//! The original network samples every path against one gas snapshot. Its
//! transported parcels then feed their explicitly bound native interfaces;
//! every regional phase/heat proposal enters one native mass/enthalpy update
//! and one flow step. Local humidity and cooling remain regional states.

use super::{
    BTreeMap, CancelGate, ClosedGasParcel, Error, GasEnergySolve, GasTransportPaths,
    GasTransportSolve, GasTransportStep, MaterialFreeSurfaceEnthalpy3, PhaseAccount,
    ThermochemicalGasModel, plumbing, poll, transport_phase_account,
};
use crate::liquid_gas::evaporation::{
    EvaporatingLiquidGasError,
    multiple::{
        DistributedLiquidGasSolve,
        regions::{
            self as phase, RegionalLiquidGasCandidate, RegionalLiquidGasHeatSurface,
            RegionalLiquidGasReport, RegionalLiquidGasSurface,
        },
    },
};
use fs_material::evaporation::gas::chamber::heat::GasLiquidChamberHeatSolve;

/// Original transport and regional phase/heat reports on one physical interval.
/// Their internal transfers already cancel; these totals combine only defects.
#[derive(Clone, Debug)]
pub struct RegionalGasTransportReport {
    pub duration_s: f64,
    pub transport: GasTransportStep,
    pub coupled: RegionalLiquidGasReport,
    pub mass_residual_kg: f64,
    pub absolute_mass_defect_kg: f64,
    pub energy_residual_j: f64,
    pub absolute_energy_defect_j: f64,
}

/// All finite gas nodes and the native liquid publish together after final gates.
pub struct RegionalGasTransportCandidate<'s, 'm, 'catalog> {
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    parcels: BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    report: RegionalGasTransportReport,
}

impl<'s, 'm, 'catalog> RegionalGasTransportCandidate<'s, 'm, 'catalog> {
    pub const fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        &self.liquid
    }
    pub const fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub const fn report(&self) -> &RegionalGasTransportReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        BTreeMap<u64, ClosedGasParcel>,
        &'catalog ThermochemicalGasModel,
        RegionalGasTransportReport,
    ) {
        (self.liquid, self.parcels, self.catalog, self.report)
    }
}

impl PhaseAccount for RegionalLiquidGasReport {
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

/// Transport the original gas network, then evaporate/condense into local gases.
///
/// Each region ID names an actual finite parcel. Different liquid cells may
/// see different gas temperature, pressure and vapor fraction; all faces of
/// one liquid cell must retain one region. Arrivals may participate in phase
/// exchange, but cannot fund outgoing paths in the same transport snapshot.
/// The complete final gas map must pass the original sampled driving-force
/// gates, so a phase-induced pressure/concentration reversal refuses the whole
/// candidate, including all liquid and gas changes.
///
/// The source catalog, per-node transport endpoint controls, shared native
/// phase controls and original material laws remain explicit. The phase
/// controls' mass/energy allowances bound the complete interval after adding
/// absolute original transport and phase/flow defects. Latent energy and
/// transported partial enthalpy are not reapplied as external heat.
///
/// This is first-order transport/phase/flow splitting on the retained native
/// timestep. Refine the whole interval. Regions are well mixed finite volumes,
/// not resolved gas CFD; this adds no regional traction, gas pressure-work,
/// finite-wall/conduction or reaction stage. Gas volumes remain fixed external
/// control volumes. Native topology and original source refusals remain active.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_regions_with_gas_transport<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    paths: GasTransportPaths<'_>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    surfaces: &[RegionalLiquidGasSurface<'_>],
    transport_controls: GasTransportSolve,
    phase_controls: DistributedLiquidGasSolve,
    gate: &CancelGate,
) -> Result<RegionalGasTransportCandidate<'s, 'm, 'catalog>, Error> {
    transporting_regions_candidate(
        liquid,
        parcels,
        paths,
        endpoints,
        transport_controls,
        phase_controls,
        gate,
        |transported| {
            phase::advance_liquid_regions_with_gas(
                liquid,
                transported,
                catalog,
                surfaces,
                phase_controls,
                gate,
            )
        },
    )
}

/// Original gas transport followed by regional reciprocal sensible heat and
/// phase exchange, then exactly one native flow step on their common interval.
///
/// Each region's heat stage updates its original gas and liquid calorics before
/// its phase stage. Surface coefficients exclude latent/carried-vapor energy.
/// Every input remains unchanged on any regional, native or final-network
/// refusal. The split accuracy, fixed volumes and model limits are those of
/// [`advance_liquid_regions_with_gas_transport`].
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_regions_with_gas_transport_and_heat<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    paths: GasTransportPaths<'_>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    surfaces: &[RegionalLiquidGasHeatSurface<'_>],
    transport_controls: GasTransportSolve,
    phase_controls: DistributedLiquidGasSolve,
    heat_controls: GasLiquidChamberHeatSolve,
    gate: &CancelGate,
) -> Result<RegionalGasTransportCandidate<'s, 'm, 'catalog>, Error> {
    transporting_regions_candidate(
        liquid,
        parcels,
        paths,
        endpoints,
        transport_controls,
        phase_controls,
        gate,
        |transported| {
            phase::advance_liquid_regions_with_gas_and_heat(
                liquid,
                transported,
                catalog,
                surfaces,
                phase_controls,
                heat_controls,
                gate,
            )
        },
    )
}

#[allow(clippy::too_many_arguments)]
fn transporting_regions_candidate<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    paths: GasTransportPaths<'_>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    transport_controls: GasTransportSolve,
    phase_controls: DistributedLiquidGasSolve,
    gate: &CancelGate,
    coupling: impl FnOnce(
        &BTreeMap<u64, ClosedGasParcel>,
    ) -> Result<
        RegionalLiquidGasCandidate<'s, 'm, 'catalog>,
        EvaporatingLiquidGasError,
    >,
) -> Result<RegionalGasTransportCandidate<'s, 'm, 'catalog>, Error> {
    poll(gate)?;
    if [
        phase_controls.mass_tolerance_kg,
        phase_controls.energy_tolerance_j,
    ]
    .into_iter()
    .any(|value| !value.is_finite() || value < 0.0)
    {
        return Err(Error::Invalid(
            "complete regional gas/liquid allowances must be finite and nonnegative",
        ));
    }
    let duration_s = liquid.flow().scaling().timestep_s();
    let (transported, transport) = plumbing::step_candidate_with_paths(
        parcels,
        paths,
        duration_s,
        endpoints,
        transport_controls,
        || gate.is_requested(),
    )?
    .into_parts();
    poll(gate)?;
    let (next_liquid, next_parcels, catalog, coupled) = coupling(&transported)?.into_parts();
    if coupled.duration_s.to_bits() != duration_s.to_bits()
        || transport.duration_s.to_bits() != duration_s.to_bits()
        || liquid.accepted_steps().checked_add(1) != Some(next_liquid.accepted_steps())
        || next_parcels.keys().ne(parcels.keys())
    {
        return Err(Error::Invalid(
            "regional transport and phase lost their gas roster or single native interval",
        ));
    }
    for (id, original) in parcels {
        poll(gate)?;
        if next_parcels[id].volume_m3().to_bits() != original.volume_m3().to_bits() {
            return Err(Error::Invalid(
                "regional phase changed an original fixed gas volume",
            ));
        }
    }
    transport.admit_endpoint_driving_forces(&next_parcels, || gate.is_requested())?;
    let account = transport_phase_account(
        &transport,
        &coupled,
        phase_controls.mass_tolerance_kg,
        phase_controls.energy_tolerance_j,
    )?;
    poll(gate)?;
    Ok(RegionalGasTransportCandidate {
        liquid: next_liquid,
        parcels: next_parcels,
        catalog,
        report: RegionalGasTransportReport {
            duration_s,
            transport,
            coupled,
            mass_residual_kg: account.mass_residual_kg,
            absolute_mass_defect_kg: account.absolute_mass_defect_kg,
            energy_residual_j: account.energy_residual_j,
            absolute_energy_defect_j: account.absolute_energy_defect_j,
        },
    })
}
