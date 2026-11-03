//! Original gas plumbing feeding native phase, heat and flow on one clock.

use super::{
    EvaporatingLiquidGasError, LiquidGasEvaporationSurface,
    conduction::{
        ConductingLiquidGasReport, advance_liquid_surfaces_with_gas_and_radiative_conduction,
    },
    multiple::DistributedLiquidGasSolve,
};
use fs_exec::CancelGate;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat, MaterialFreeSurfaceEnthalpy3,
    conduction::{
        FreeSurfaceBodyContact, FreeSurfaceBodyInterface, FreeSurfaceConductionControls,
        FreeSurfaceConductionFace, FreeSurfaceGasFace, FreeSurfaceRadiationFace,
        FreeSurfaceThermalBody,
    },
};
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergySolve,
        parcel::{
            ClosedGasParcel,
            transfer::species::network::transport::{
                self as plumbing, GasTransportError, GasTransportPaths, GasTransportSolve,
                GasTransportStep,
            },
        },
    },
};
use std::collections::BTreeMap;

#[derive(Debug)]
pub enum TransportingLiquidGasError {
    Transport(GasTransportError),
    Coupled(EvaporatingLiquidGasError),
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
impl core::fmt::Display for TransportingLiquidGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "gas plumbing/native liquid phase and heat: {self:?}")
    }
}
impl std::error::Error for TransportingLiquidGasError {}
impl From<GasTransportError> for TransportingLiquidGasError {
    fn from(value: GasTransportError) -> Self {
        match value {
            GasTransportError::Cancelled => Self::Cancelled,
            value => Self::Transport(value),
        }
    }
}
impl From<EvaporatingLiquidGasError> for TransportingLiquidGasError {
    fn from(value: EvaporatingLiquidGasError) -> Self {
        match value {
            EvaporatingLiquidGasError::Cancelled => Self::Cancelled,
            value => Self::Coupled(value),
        }
    }
}

/// Both original stage reports and the complete source-relative budget.
/// Transport entropy is local to its closed gas stage; it is not a claim about
/// the subsequent open gas/liquid thermal subsystem.
#[derive(Clone, Debug)]
pub struct TransportingLiquidGasReport {
    pub duration_s: f64,
    pub atmosphere_id: u64,
    pub transport: GasTransportStep,
    pub coupled: ConductingLiquidGasReport,
    pub mass_residual_kg: f64,
    pub absolute_mass_defect_kg: f64,
    pub energy_residual_j: f64,
    pub absolute_energy_defect_j: f64,
}

/// Every finite gas node, native liquid and finite wall publishes together.
pub struct TransportingLiquidGasCandidate<'s, 'm, 'catalog, 'b> {
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    parcels: BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    bodies: Vec<FreeSurfaceThermalBody<'b>>,
    report: TransportingLiquidGasReport,
}
impl<'s, 'm, 'catalog, 'b> TransportingLiquidGasCandidate<'s, 'm, 'catalog, 'b> {
    pub fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        &self.liquid
    }
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub fn gas(&self) -> &ClosedGasParcel {
        &self.parcels[&self.report.atmosphere_id]
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub fn bodies(&self) -> &[FreeSurfaceThermalBody<'b>] {
        &self.bodies
    }
    pub const fn report(&self) -> &TransportingLiquidGasReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        BTreeMap<u64, ClosedGasParcel>,
        &'catalog ThermochemicalGasModel,
        Vec<FreeSurfaceThermalBody<'b>>,
        TransportingLiquidGasReport,
    ) {
        (
            self.liquid,
            self.parcels,
            self.catalog,
            self.bodies,
            self.report,
        )
    }
}

type Error = TransportingLiquidGasError;
fn poll(gate: &CancelGate) -> Result<(), Error> {
    if gate.is_requested() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

/// Feed one native liquid atmosphere through the original finite-gas network,
/// then exchange phase/heat and advance native flow ONCE over the same dt.
///
/// All plumbing paths sample the unchanged input gas map and share its original
/// gross withdrawal allowance. The nominated atmosphere's exact transported
/// parcel then enters the original phase/conduction owner. Its catalog must
/// cover the original source identities of the transported composition; no
/// replacement species, transport law or caloric reference is inferred.
///
/// Finite wall contact, gas films, prescribed liquid/wall radiation, attached
/// electrical heat and optional physical viscous heat remain with the existing
/// native thermal owner. Optional momentum loading uses the post-thermal gas.
/// The actual final COMPLETE gas map must still satisfy the original sampled
/// transport driving-force and source gates. A phase/thermal pressure reversal
/// therefore rejects all candidate nodes, liquid and walls together.
///
/// This is explicit first-order transport/phase/conduction/flow splitting.
/// Arrivals may participate in the subsequent phase stage, but cannot fund an
/// outgoing plumbing path in the same transport snapshot. It is not a joint
/// implicit gas-flow/phase solve, spatial gas CFD, jet momentum or volume-work
/// model. Refine the complete interval for splitting and transport accuracy.
///
/// The phase controls' mass/energy allowances also bound the complete interval:
/// original absolute transport defects are added to original coupled defects.
/// Internal transported partial enthalpy and phase latent heat are never added
/// again as external heat. Gas-only transport entropy remains stage-local.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_surfaces_with_gas_transport<'s, 'm, 'catalog, 'b>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    atmosphere_id: u64,
    catalog: &'catalog ThermochemicalGasModel,
    paths: GasTransportPaths<'_>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    surfaces: &[LiquidGasEvaporationSurface<'_>],
    bodies: &[FreeSurfaceThermalBody<'b>],
    faces: &[FreeSurfaceConductionFace],
    heat: &[FreeSurfaceCellHeat],
    contacts: &[FreeSurfaceBodyContact],
    interfaces: &[FreeSurfaceBodyInterface<'_>],
    exposed: &[FreeSurfaceGasFace],
    radiation: &[FreeSurfaceRadiationFace],
    transport_controls: GasTransportSolve,
    phase_controls: DistributedLiquidGasSolve,
    thermal_controls: FreeSurfaceConductionControls,
    maximum_density_excursion: Option<f64>,
    with_viscous_heating: bool,
    gate: &CancelGate,
) -> Result<TransportingLiquidGasCandidate<'s, 'm, 'catalog, 'b>, Error> {
    poll(gate)?;
    if !parcels.contains_key(&atmosphere_id) {
        return Err(Error::Invalid(
            "native atmosphere must name an original finite gas node",
        ));
    }
    let gas_controls = *endpoints.get(&atmosphere_id).ok_or(Error::Invalid(
        "native atmosphere requires its original gas endpoint controls",
    ))?;
    if [
        phase_controls.mass_tolerance_kg,
        phase_controls.energy_tolerance_j,
    ]
    .into_iter()
    .any(|value| !value.is_finite() || value < 0.0)
    {
        return Err(Error::Invalid(
            "complete gas/liquid mass and energy allowances must be finite and nonnegative",
        ));
    }
    let duration_s = liquid.flow().scaling().timestep_s();
    let (mut next_parcels, transport) = plumbing::step_candidate_with_paths(
        parcels,
        paths,
        duration_s,
        endpoints,
        transport_controls,
        || gate.is_requested(),
    )?
    .into_parts();
    let transported = next_parcels.get(&atmosphere_id).ok_or(Error::Invalid(
        "original gas transport lost the native atmosphere node",
    ))?;
    let candidate = advance_liquid_surfaces_with_gas_and_radiative_conduction(
        liquid,
        transported,
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
        gate,
    )?;
    let (next_liquid, next_gas, catalog, next_bodies, coupled) = candidate.into_parts();
    if coupled.phase.duration_s.to_bits() != duration_s.to_bits()
        || transport.duration_s.to_bits() != duration_s.to_bits()
        || liquid.accepted_steps().checked_add(1) != Some(next_liquid.accepted_steps())
        || next_gas.volume_m3().to_bits() != parcels[&atmosphere_id].volume_m3().to_bits()
    {
        return Err(Error::Invalid(
            "gas plumbing and native liquid lost their one interval or fixed atmosphere volume",
        ));
    }
    next_parcels.insert(atmosphere_id, next_gas);
    transport.admit_endpoint_driving_forces(&next_parcels, || gate.is_requested())?;

    let original = &transport.exchange;
    let transport_mass_defect = original
        .parcels
        .iter()
        .map(|row| {
            row.mass_residual_kg
                .abs()
                .max(row.maximum_species_mass_residual_kg)
        })
        .sum::<f64>()
        .max(original.total_mass_residual_kg.abs())
        .max(original.maximum_species_mass_residual_kg);
    let transport_energy_defect = original
        .parcels
        .iter()
        .map(|row| row.energy_residual_j.abs())
        .sum::<f64>()
        .max(original.total_energy_residual_j.abs());
    let mass_residual_kg = original.total_mass_residual_kg + coupled.phase.mass_residual_kg;
    let absolute_mass_defect_kg = transport_mass_defect + coupled.phase.absolute_mass_defect_kg;
    let energy_residual_j = original.total_energy_residual_j + coupled.energy_residual_j;
    let absolute_energy_defect_j = transport_energy_defect + coupled.absolute_energy_defect_j;
    if [
        mass_residual_kg,
        absolute_mass_defect_kg,
        energy_residual_j,
        absolute_energy_defect_j,
    ]
    .into_iter()
    .any(|value| !value.is_finite())
    {
        return Err(Error::Invalid(
            "complete gas plumbing/phase/thermal account is unrepresentable",
        ));
    }
    if absolute_mass_defect_kg > phase_controls.mass_tolerance_kg {
        return Err(Error::MassBalance {
            absolute_defect_kg: absolute_mass_defect_kg,
            tolerance_kg: phase_controls.mass_tolerance_kg,
        });
    }
    if absolute_energy_defect_j > phase_controls.energy_tolerance_j {
        return Err(Error::EnergyBalance {
            absolute_defect_j: absolute_energy_defect_j,
            tolerance_j: phase_controls.energy_tolerance_j,
        });
    }
    poll(gate)?;
    Ok(TransportingLiquidGasCandidate {
        liquid: next_liquid,
        parcels: next_parcels,
        catalog,
        bodies: next_bodies,
        report: TransportingLiquidGasReport {
            duration_s,
            atmosphere_id,
            transport,
            coupled,
            mass_residual_kg,
            absolute_mass_defect_kg,
            energy_residual_j,
            absolute_energy_defect_j,
        },
    })
}
