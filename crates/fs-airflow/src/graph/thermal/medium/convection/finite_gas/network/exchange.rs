//! Mixed material gas heat laws on the original reacting/radiative assembly.
//! Each patch selects its actual law; no material-name or phase-name dispatch.
use super::*;
use super::super::NaturalGasSurface;
use fs_couple::cavity_phs::piston::thermal::ConductanceSample;
use fs_couple::cavity_phs::piston::thermal::boundary::network::material::GasConductionGeometry;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::{
    MaterialContact, exchange::{MaterialExchangeSample, PreparedMaterialExchangeNetwork, gas_law::GasSurface},
    reactive::{ChamberReaction, ReactiveExchangeCandidate, ReactiveExchangeError},
    state::{ChamberDrive, MaterialExchangeState},
    state::flow::{MaterialFlowAdvance, MaterialGasLawFlowError},
};
use fs_exec::Cx;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::{
    GasNetworkLink, GasNetworkSolve, GasNetworkStep,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::SnapshotNozzleLink,
    species::{diffusion::network::DiffusionNetworkLink,
        network::transport::{GasTransportPaths, GasTransportSolve, GasTransportStep}},
};

/// Constitutive choice for a disjoint physical patch, not an object preset.
#[derive(Clone, Copy, Debug)]
pub enum GasHeatLaw {
    /// Original source-backed gas and body layers in series.
    Conduction(GasConductionGeometry),
    /// Original forced correlation at the actual gas energy midpoint.
    ForcedConvection(ForcedGasSurface),
    /// Buoyant vertical-plate transfer from the same material gas and wall state.
    NaturalConvection(NaturalGasSurface),
}
/// One patch using the same chamber/cell IDs as the complete assembly.
#[derive(Clone, Copy, Debug)]
pub struct ExchangeContact {
    pub contact: Contact,
    pub law: GasHeatLaw,
}
/// Preserve the actual source result; no synthetic common material receipt.
#[derive(Debug)]
pub enum GasHeatSample {
    Conduction(ResolvedMaterialStatePoint),
    ForcedConvection(FiniteGasConvectionSample),
    NaturalConvection(FiniteGasConvectionSample),
}
/// Preserve the original domain refusal, not an invented substitute coefficient.
#[derive(Debug)]
pub enum GasHeatError {
    Conduction(MaterialNetworkError),
    Convection(MaterialConvectionError),
}
impl core::fmt::Display for GasHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material gas heat law: {self:?}")
    }
}
impl std::error::Error for GasHeatError {}

impl GasHeatLaw {
    fn area(self) -> Result<f64, &'static str> {
        match self {
            Self::Conduction(geometry) => {
                if ![geometry.area_m2, geometry.gas_distance_m, geometry.cell_distance_m]
                    .iter().all(|v| v.is_finite())
                    || geometry.area_m2 < 0.0 || geometry.gas_distance_m <= 0.0
                    || geometry.cell_distance_m <= 0.0 {
                    return Err("valid conductive patch area and positive layer lengths required");
                }
                Ok(geometry.area_m2)
            }
            Self::ForcedConvection(surface) => {
                surface.validate().map_err(|_| "invalid forced-convection geometry")?;
                Ok(surface.area_m2)
            }
            Self::NaturalConvection(surface) => {
                surface.validate().map_err(|_| "invalid natural-convection geometry")?;
                Ok(surface.area_m2)
            }
        }
    }

    fn resolve(
        self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, contact: Contact,
        thermal: &DepositionSimulation<'_>, source: &MaterialExchangeSample,
        gas: &ResolvedThermochemicalGasState,
    ) -> Result<ConductanceSample<GasHeatSample>, GasHeatError> {
        if let Self::Conduction(geometry) = self {
            let sample = prepared.resolve_conductive_gas_contact(
                MaterialContact { contact, geometry }, thermal, source, gas)
                .map_err(GasHeatError::Conduction)?;
            return Ok(ConductanceSample { conductance_w_k: sample.conductance_w_k,
                source: GasHeatSample::Conduction(sample.source) });
        }
        let phase = source.bulk.cells.iter().find(|c| c.cell == contact.cell)
            .ok_or_else(|| GasHeatError::Convection(
                invalid("gas patch lacks its actual cell phase")))?.phase;
        let (sample, natural) = match self {
            Self::ForcedConvection(surface) => (surface.resolve(gas, phase), false),
            Self::NaturalConvection(surface) => (surface.resolve(gas, phase), true),
            Self::Conduction(_) => unreachable!("conduction returned above"),
        };
        let sample = sample.map_err(GasHeatError::Convection)?;
        Ok(ConductanceSample { conductance_w_k: sample.conductance_w_k,
            source: if natural { GasHeatSample::NaturalConvection(sample.source) }
                else { GasHeatSample::ForcedConvection(sample.source) } })
    }
}

fn flow_patches<'a>(
    cx: &Cx<'_>, patches: &'a [ExchangeContact], controls: NetworkSolve,
) -> Result<(BTreeMap<u64, &'a ExchangeContact>, Vec<GasSurface>), MaterialGasLawFlowError<GasHeatError>> {
    if patches.len() > controls.maximum_contacts {
        return Err(MaterialGasLawFlowError::Invalid("mixed heat-law contact capacity exceeded"));
    }
    let mut mapped = BTreeMap::new();
    let mut surfaces = Vec::with_capacity(patches.len());
    for patch in patches {
        cx.checkpoint().map_err(|_| MaterialGasLawFlowError::Cancelled)?;
        let area = patch.law.area().map_err(MaterialGasLawFlowError::Invalid)?;
        if mapped.insert(patch.contact.id, patch).is_some() {
            return Err(MaterialGasLawFlowError::Invalid("duplicate mixed heat-law contact ID"));
        }
        surfaces.push(GasSurface { contact: patch.contact, area_m2: area });
    }
    Ok((mapped, surfaces))
}

/// Compose mixed gas heat laws with the SAME declared reactions, electrical
/// contact losses, thermal interface resistance, radiation, finite cell storage
/// and wall motion. Empty reactions retain the original nonreacting physics.
///
/// Forced and natural patches consume the same actual source gas and caloric
/// wall endpoint. Natural convection needs no authored flow speed: the shared
/// vertical-plate law resolves buoyant transfer from those temperatures and
/// the gas's ideal EOS, viscosity, conductivity and heat capacity. Its bulk-gas
/// property convention and unconfined vertical-plate assumptions must fit the
/// physical patch; a chamber ID does not establish enclosure applicability.
///
/// Each chemical stage runs once on its original gas (including optional
/// chemical subcycling). Then the existing coupled owner solves every thermal
/// and mechanical trial against that reacted state. Each convection patch sees
/// the actual trial gas; each conduction patch uses the complete pinned cell
/// bundle. No rate is frozen outside its chemical owner, repeated in Newton,
/// converted into extra heat, or resampled by this adapter.
///
/// All conduction AND convection patch areas share one exposure check, across
/// all chambers. Radiation may share a physical patch; two heat laws may not
/// silently claim overlapping patches. Non-overlap is a caller geometry claim.
/// Conductivity is required only where an actual conduction law needs it.
///
/// The returned candidate owns circuit, cells, gases and momenta together.
/// Preserve its original source records and publish all parts at once. This
/// adapter does not add gas mass transport, pump work, wall deformation, a
/// boundary-layer solve or a stronger temporal order than the original assembly.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m>(
    prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>,
    thermal: &DepositionSimulation<'m>, target_s: f64,
    circuit: &MaterialCircuit, ports: &[f64], chambers: &[Chamber<'_>],
    reactions: &[ChamberReaction<'_>], patches: &[ExchangeContact],
    controls: NetworkSolve, guess_j: &[f64],
) -> Result<ReactiveExchangeCandidate<'m, GasHeatSample>, ReactiveExchangeError<GasHeatError>> {
    let check = || -> Result<(), ReactiveExchangeError<GasHeatError>> {
        cx.checkpoint().map_err(|_| ReactiveExchangeError::Cancelled)
    };
    check()?;
    if patches.len() > controls.maximum_contacts || chambers.len() > controls.maximum_chambers {
        return Err(ReactiveExchangeError::Invalid("mixed heat-law chamber/contact capacity exceeded"));
    }
    let mut mapped = BTreeMap::new();
    let mut surfaces = Vec::with_capacity(patches.len());
    for patch in patches {
        check()?;
        let area = patch.law.area().map_err(ReactiveExchangeError::Invalid)?;
        if mapped.insert(patch.contact.id, patch).is_some() {
            return Err(ReactiveExchangeError::Invalid("duplicate mixed heat-law contact ID"));
        }
        surfaces.push(GasSurface { contact: patch.contact, area_m2: area });
    }
    prepared.reactive_with_gas_law_candidate(cx, thermal, target_s, circuit,
        ports, chambers, reactions, &surfaces, controls, guess_j,
        |contact, state, source, gas| mapped[&contact.id].law.resolve(
            prepared, contact, state, source, gas))
}


/// Preview conservative nozzle transport followed by the same mixed conductive
/// and convective heat laws, optional chemistry, finite electrical/thermal
/// storage and wall motion owned by a live MaterialExchangeState.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn flow_step_candidate<'m>(
    state: &MaterialExchangeState<'m>,
    prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>,
    target_s: f64, ports: &[f64], drives: &[ChamberDrive],
    reactions: &[ChamberReaction<'_>], patches: &[ExchangeContact],
    links: &[GasNetworkLink], flow_controls: GasNetworkSolve,
    controls: NetworkSolve, guess_j: &[f64],
) -> Result<(MaterialExchangeState<'m>, MaterialFlowAdvance<GasNetworkStep, GasHeatSample>),
    MaterialGasLawFlowError<GasHeatError>> {
    let (mapped, surfaces) = flow_patches(cx, patches, controls)?;
    state.flow_step_with_gas_law_candidate(prepared,cx,target_s,ports,drives,
        reactions,&surfaces,links,flow_controls,controls,guess_j,
        |contact, thermal, source, gas| mapped[&contact.id].law.resolve(
            prepared, contact, thermal, source, gas))
}

/// Publish the complete transport/reaction/mixed-heat assembly atomically.
#[allow(clippy::too_many_arguments)]
pub fn advance_with_flow<'m>(
    state: &mut MaterialExchangeState<'m>,
    prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>,
    target_s: f64, ports: &[f64], drives: &[ChamberDrive],
    reactions: &[ChamberReaction<'_>], patches: &[ExchangeContact],
    links: &[GasNetworkLink], flow_controls: GasNetworkSolve,
    controls: NetworkSolve, guess_j: &[f64],
) -> Result<MaterialFlowAdvance<GasNetworkStep, GasHeatSample>,
    MaterialGasLawFlowError<GasHeatError>> {
    let (next,report)=flow_step_candidate(state,prepared,cx,target_s,ports,drives,
        reactions,patches,links,flow_controls,controls,guess_j)?;
    cx.checkpoint().map_err(|_|MaterialGasLawFlowError::Cancelled)?;
    *state=next;
    Ok(report)
}

/// Preview same-snapshot pressure flow AND molecular mixing, then the original
/// reacting/electrical/thermal/mechanical assembly with mixed gas heat laws.
///
/// Both transport mechanisms consume one original species budget, reconstruct
/// each chamber once and carry partial enthalpy once. Conduction, forced and
/// natural patches all receive the actual transported/reacted trial gas. Thus
/// diffusion changes heat transfer, storage and sound speed through shared
/// material data even at equal total pressure, with no ambient composition cache.
///
/// Nozzles and diffusion faces describe distinct physical paths. Empty lists
/// disable that mechanism. The existing explicit first-order split, diffusion
/// temperature/pressure domains and convection geometry restrictions apply;
/// forced speeds are still prescribed, not inferred from nozzle jet velocity.
/// Every physical endpoint and its port history is returned together or none.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn transport_step_candidate<'m>(
    state: &MaterialExchangeState<'m>,
    prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>,
    target_s: f64, ports: &[f64], drives: &[ChamberDrive],
    reactions: &[ChamberReaction<'_>], patches: &[ExchangeContact],
    nozzles: &[SnapshotNozzleLink], diffusion: &[DiffusionNetworkLink<'_>],
    transport_controls: GasTransportSolve, controls: NetworkSolve, guess_j: &[f64],
) -> Result<(MaterialExchangeState<'m>, MaterialFlowAdvance<GasTransportStep, GasHeatSample>),
    MaterialGasLawFlowError<GasHeatError>> {
    transport_paths_step_candidate(state, prepared, cx, target_s, ports, drives,
        reactions, patches, GasTransportPaths {
            nozzles, diffusion, ..GasTransportPaths::default()
        }, transport_controls, controls, guess_j)
}

/// Preview nozzle, porous, molecular and membrane transport before the same
/// reacting assembly and mixed heat laws. Every path spends one shared original
/// species allowance; changed composition supplies the actual thermal trial gas.
///
/// Darcy and membrane constitutive domains remain explicit. Paths describe
/// separately authored exposed areas. Neither porous drag nor carried species
/// enthalpy is added again as wall heat. Final driving-force admission encloses
/// the thermal/mechanical step.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn transport_paths_step_candidate<'m>(
    state: &MaterialExchangeState<'m>,
    prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>,
    target_s: f64, ports: &[f64], drives: &[ChamberDrive],
    reactions: &[ChamberReaction<'_>], patches: &[ExchangeContact],
    paths: GasTransportPaths<'_>, transport_controls: GasTransportSolve,
    controls: NetworkSolve, guess_j: &[f64],
) -> Result<(MaterialExchangeState<'m>, MaterialFlowAdvance<GasTransportStep, GasHeatSample>),
    MaterialGasLawFlowError<GasHeatError>> {
    let (mapped, surfaces) = flow_patches(cx, patches, controls)?;
    state.transport_paths_step_with_gas_law_candidate(prepared, cx, target_s, ports,
        drives, reactions, &surfaces, paths, transport_controls,
        controls, guess_j, |contact, thermal, source, gas|
            mapped[&contact.id].law.resolve(prepared, contact, thermal, source, gas))
}

/// Publish combined species transport, mixed gas heat laws and every original
/// coupled owner atomically. A late material/domain/solve refusal changes none.
#[allow(clippy::too_many_arguments)]
pub fn advance_with_transport<'m>(
    state: &mut MaterialExchangeState<'m>,
    prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>,
    target_s: f64, ports: &[f64], drives: &[ChamberDrive],
    reactions: &[ChamberReaction<'_>], patches: &[ExchangeContact],
    nozzles: &[SnapshotNozzleLink], diffusion: &[DiffusionNetworkLink<'_>],
    transport_controls: GasTransportSolve, controls: NetworkSolve, guess_j: &[f64],
) -> Result<MaterialFlowAdvance<GasTransportStep, GasHeatSample>,
    MaterialGasLawFlowError<GasHeatError>> {
    let (next, report) = transport_step_candidate(state, prepared, cx, target_s,
        ports, drives, reactions, patches, nozzles, diffusion, transport_controls,
        controls, guess_j)?;
    cx.checkpoint().map_err(|_| MaterialGasLawFlowError::Cancelled)?;
    *state = next;
    Ok(report)
}

/// Publish all gas transport and the complete reacting/thermal/mechanical owner
/// together. Any late refusal leaves every original state and history unchanged.
#[allow(clippy::too_many_arguments)]
pub fn advance_with_transport_paths<'m>(
    state: &mut MaterialExchangeState<'m>,
    prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>,
    target_s: f64, ports: &[f64], drives: &[ChamberDrive],
    reactions: &[ChamberReaction<'_>], patches: &[ExchangeContact],
    paths: GasTransportPaths<'_>, transport_controls: GasTransportSolve,
    controls: NetworkSolve, guess_j: &[f64],
) -> Result<MaterialFlowAdvance<GasTransportStep, GasHeatSample>,
    MaterialGasLawFlowError<GasHeatError>> {
    let (next, report) = transport_paths_step_candidate(state, prepared, cx,
        target_s, ports, drives, reactions, patches, paths, transport_controls,
        controls, guess_j)?;
    cx.checkpoint().map_err(|_| MaterialGasLawFlowError::Cancelled)?;
    *state = next;
    Ok(report)
}
