//! Mixed material gas heat laws on the original reacting/radiative assembly.
//! Each patch selects its actual law; no material-name or phase-name dispatch.
use super::*;
use fs_couple::cavity_phs::piston::thermal::ConductanceSample;
use fs_couple::cavity_phs::piston::thermal::boundary::network::material::GasConductionGeometry;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::{
    MaterialContact, exchange::{PreparedMaterialExchangeNetwork, gas_law::GasSurface},
    reactive::{ChamberReaction, ReactiveExchangeCandidate, ReactiveExchangeError},
};
use fs_exec::Cx;
use fs_material::state_point::ResolvedMaterialStatePoint;

/// Constitutive choice for a disjoint physical patch, not an object preset.
#[derive(Clone, Copy, Debug)]
pub enum GasHeatLaw {
    /// Original source-backed gas and body layers in series.
    Conduction(GasConductionGeometry),
    /// Original forced correlation at the actual gas energy midpoint.
    ForcedConvection(ForcedGasSurface),
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

/// Compose mixed gas heat laws with the SAME declared reactions, electrical
/// contact losses, thermal interface resistance, radiation, finite cell storage
/// and wall motion. Empty reactions retain the original nonreacting physics.
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
        let area = match patch.law {
            GasHeatLaw::Conduction(geometry) => {
                if ![geometry.gas_distance_m, geometry.cell_distance_m].iter()
                    .all(|v| v.is_finite() && *v > 0.0) {
                    return Err(ReactiveExchangeError::Invalid("positive finite conductive layer lengths required"));
                }
                geometry.area_m2
            }
            GasHeatLaw::ForcedConvection(surface) => {
                surface.validate().map_err(|_| ReactiveExchangeError::Invalid("invalid forced-convection geometry"))?;
                surface.area_m2
            }
        };
        if mapped.insert(patch.contact.id, patch).is_some() {
            return Err(ReactiveExchangeError::Invalid("duplicate mixed heat-law contact ID"));
        }
        surfaces.push(GasSurface { contact: patch.contact, area_m2: area });
    }
    prepared.reactive_with_gas_law_candidate(cx, thermal, target_s, circuit,
        ports, chambers, reactions, &surfaces, controls, guess_j,
        |contact, state, source, gas| {
            match mapped[&contact.id].law {
                GasHeatLaw::Conduction(geometry) => {
                    let sample = prepared.resolve_conductive_gas_contact(
                        MaterialContact { contact, geometry }, state, source, gas)
                        .map_err(GasHeatError::Conduction)?;
                    Ok(ConductanceSample { conductance_w_k: sample.conductance_w_k,
                        source: GasHeatSample::Conduction(sample.source) })
                }
                GasHeatLaw::ForcedConvection(surface) => {
                    let phase = source.bulk.cells.iter().find(|c| c.cell == contact.cell)
                        .ok_or_else(|| GasHeatError::Convection(invalid("gas patch lacks its actual cell phase")))?.phase;
                    let sample = surface.resolve(gas, phase).map_err(GasHeatError::Convection)?;
                    Ok(ConductanceSample { conductance_w_k: sample.conductance_w_k,
                        source: GasHeatSample::ForcedConvection(sample.source) })
                }
            }
        })
}
