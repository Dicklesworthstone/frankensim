//! Original gas-network transport feeding a collective spherical atmosphere.
//!
//! One original transport proposal advances every declared gas path. Every
//! droplet film then samples the ACTUAL transported gas in its explicitly
//! named region, before one original shared phase/optional sensible-heat solve.
//! The accepted gas returns to the same complete network map, and the original
//! transport owner rechecks its driving forces at the final coupled endpoint.
//! A late refusal publishes neither the transported gases nor any droplet.
//!
//! The common interval is first-order transport then heat/phase splitting;
//! source and step-size refusals require retry from the original states. Gas
//! volumes remain prescribed and fixed. The caller supplies droplet relative
//! speeds, heat coefficients and interval refinement; gas transport does not
//! infer droplet trajectories, drag, entrainment, shielding or moving-volume
//! work. The transport entropy receipt describes that substep only.

/// Adaptive physical time for the complete original gas and droplet network.
pub mod adaptive;

use std::collections::BTreeMap;

use super::{
    ClosedGasParcel, EvaporationError, GasFilmError, SourceDropletFilmStep,
    SourceGasDropletChamberCandidate, SourceGasDropletHeatSurface, SourceGasDropletSurface,
    SphericalDropletFilmError, ThermochemicalGasModel,
    advance_droplets_with_gas as advance_chamber,
    advance_droplets_with_gas_and_heat as advance_chamber_with_heat,
};
use fs_material::evaporation::gas::{
    GasEvaporationError,
    chamber::{GasLiquidChamberSolve, heat::GasLiquidChamberHeatSolve},
    droplet::chamber::{GasDropletChamberReport, GasDropletSurfaceStep},
};
use fs_material::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::transfer::species::network::transport::{
        self as plumbing, GasTransportError, GasTransportPaths, GasTransportSolve, GasTransportStep,
    },
};

/// Original per-stage controls plus one extensive whole-interval allowance.
/// The phase owner retains its explicit gas caloric bracket; transport uses
/// the supplied per-node endpoint controls. Neither stage enlarges the other.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasDropletTransportSolve {
    pub transport: GasTransportSolve,
    pub phase: GasLiquidChamberSolve,
    pub mass_tolerance_kg: f64,
    pub energy_tolerance_j: f64,
}

/// Original source/transport failures remain identifiable, including retries.
#[derive(Debug)]
pub enum DropletGasTransportError {
    Invalid(&'static str),
    Transport(GasTransportError),
    Droplets(SphericalDropletFilmError),
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
impl core::fmt::Display for DropletGasTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "gas transport with collective spherical droplets: {self:?}"
        )
    }
}
impl std::error::Error for DropletGasTransportError {}
impl From<GasTransportError> for DropletGasTransportError {
    fn from(error: GasTransportError) -> Self {
        match error {
            GasTransportError::Cancelled => Self::Cancelled,
            other => Self::Transport(other),
        }
    }
}
impl From<SphericalDropletFilmError> for DropletGasTransportError {
    fn from(error: SphericalDropletFilmError) -> Self {
        match error {
            SphericalDropletFilmError::Film(GasFilmError::Cancelled)
            | SphericalDropletFilmError::Liquid(EvaporationError::Cancelled)
            | SphericalDropletFilmError::Phase(GasEvaporationError::Cancelled) => Self::Cancelled,
            other => Self::Droplets(other),
        }
    }
}
use DropletGasTransportError as Error;

/// Original transport and droplet receipts share one physical interval.
/// Carried species enthalpy, sensible heat and latent subaccounts are already
/// applied by their original owners; only defects are combined here.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletTransportReport {
    pub duration_s: f64,
    pub gas_region: u64,
    pub transport: GasTransportStep,
    pub coupled: GasDropletChamberReport,
    pub mass_residual_kg: f64,
    pub absolute_mass_defect_kg: f64,
    pub energy_residual_j: f64,
    pub absolute_energy_defect_j: f64,
}

/// All physical inventories and source receipts publish atomically.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletTransportCandidate<'catalog> {
    droplets: Vec<GasDropletSurfaceStep>,
    parcels: BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    films: Vec<SourceDropletFilmStep>,
    report: GasDropletTransportReport,
}
impl<'catalog> GasDropletTransportCandidate<'catalog> {
    pub fn droplets(&self) -> &[GasDropletSurfaceStep] {
        &self.droplets
    }
    pub const fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub fn films(&self) -> &[SourceDropletFilmStep] {
        &self.films
    }
    pub const fn report(&self) -> &GasDropletTransportReport {
        &self.report
    }
    #[allow(clippy::type_complexity)]
    pub fn into_parts(
        self,
    ) -> (
        Vec<GasDropletSurfaceStep>,
        BTreeMap<u64, ClosedGasParcel>,
        &'catalog ThermochemicalGasModel,
        Vec<SourceDropletFilmStep>,
        GasDropletTransportReport,
    ) {
        (
            self.droplets,
            self.parcels,
            self.catalog,
            self.films,
            self.report,
        )
    }
}

/// Advance original gas plumbing, then all sourced droplets in one named gas.
/// Empty paths reduce to the source chamber while retaining the complete gas
/// roster. Empty droplet rosters still run the original gas transport owner.
#[allow(clippy::too_many_arguments)]
pub fn advance_droplets_with_gas_transport<'catalog>(
    surfaces: &[SourceGasDropletSurface<'_, '_>],
    gas_region: u64,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    paths: GasTransportPaths<'_>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    duration_s: f64,
    controls: GasDropletTransportSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasDropletTransportCandidate<'catalog>, Error> {
    advance(
        gas_region,
        parcels,
        paths,
        endpoints,
        duration_s,
        controls,
        |gas, check| advance_chamber(surfaces, gas, catalog, duration_s, controls.phase, check),
        cancelled,
    )
}

/// The same transaction with one collective reciprocal sensible-heat solve.
/// Sourced geometry and transport freeze after gas transport and before heat;
/// the original chamber then owns joint heat followed by joint phase transfer.
#[allow(clippy::too_many_arguments)]
pub fn advance_droplets_with_gas_transport_and_heat<'catalog>(
    surfaces: &[SourceGasDropletHeatSurface<'_, '_>],
    gas_region: u64,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    paths: GasTransportPaths<'_>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    duration_s: f64,
    controls: GasDropletTransportSolve,
    heat_controls: GasLiquidChamberHeatSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasDropletTransportCandidate<'catalog>, Error> {
    advance(
        gas_region,
        parcels,
        paths,
        endpoints,
        duration_s,
        controls,
        |gas, check| {
            advance_chamber_with_heat(
                surfaces,
                gas,
                catalog,
                duration_s,
                controls.phase,
                heat_controls,
                check,
            )
        },
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance<'catalog>(
    gas_region: u64,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    paths: GasTransportPaths<'_>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    duration_s: f64,
    controls: GasDropletTransportSolve,
    coupling: impl FnOnce(
        &ClosedGasParcel,
        &mut dyn FnMut() -> bool,
    ) -> Result<
        SourceGasDropletChamberCandidate<'catalog>,
        SphericalDropletFilmError,
    >,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasDropletTransportCandidate<'catalog>, Error> {
    poll(&mut cancelled)?;
    let original_gas = parcels.get(&gas_region).ok_or(Error::Invalid(
        "spherical droplets must name an original finite gas region",
    ))?;
    if [controls.mass_tolerance_kg, controls.energy_tolerance_j]
        .into_iter()
        .any(|value| !value.is_finite() || value < 0.0)
    {
        return Err(Error::Invalid(
            "whole gas/droplet mass and energy allowances must be finite and nonnegative",
        ));
    }
    let (mut next_parcels, transport) = plumbing::step_candidate_with_paths(
        parcels,
        paths,
        duration_s,
        endpoints,
        controls.transport,
        &mut cancelled,
    )?
    .into_parts();
    let transported = next_parcels.get(&gas_region).ok_or(Error::Invalid(
        "original gas transport lost the droplet atmosphere",
    ))?;
    let (original, films) = coupling(transported, &mut cancelled)?.into_parts();
    let (droplets, gas, catalog, coupled) = original.into_parts();
    if transport.duration_s.to_bits() != duration_s.to_bits()
        || coupled.exchange.duration_s.to_bits() != duration_s.to_bits()
        || coupled.phase.duration_s.to_bits() != duration_s.to_bits()
        || gas.volume_m3().to_bits() != original_gas.volume_m3().to_bits()
    {
        return Err(Error::Invalid(
            "gas transport and droplets must retain one interval and fixed atmosphere volume",
        ));
    }
    // Move the sole accepted chamber gas into the already-owned network map.
    // This neither reconstructs nor clones another physical atmosphere.
    next_parcels.insert(gas_region, gas);
    transport.admit_endpoint_driving_forces(&next_parcels, &mut cancelled)?;

    let mut local_mass = 0.0;
    let mut local_energy = 0.0;
    for row in &transport.exchange.parcels {
        poll(&mut cancelled)?;
        local_mass += row
            .mass_residual_kg
            .abs()
            .max(row.maximum_species_mass_residual_kg);
        local_energy += row.energy_residual_j.abs();
    }
    let transport_mass_defect = local_mass
        .max(transport.exchange.total_mass_residual_kg.abs())
        .max(transport.exchange.maximum_species_mass_residual_kg);
    let transport_energy_defect =
        local_energy.max(transport.exchange.total_energy_residual_j.abs());
    let mass_residual_kg =
        transport.exchange.total_mass_residual_kg + coupled.exchange.mass_residual_kg;
    let absolute_mass_defect_kg =
        (transport_mass_defect + coupled.sum_absolute_mass_residual_kg).max(mass_residual_kg.abs());
    let energy_residual_j =
        transport.exchange.total_energy_residual_j + coupled.exchange.energy_residual_j;
    let absolute_energy_defect_j = (transport_energy_defect
        + coupled.exchange.sum_absolute_energy_residual_j)
        .max(energy_residual_j.abs());
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
            "complete gas transport/droplet account is unrepresentable",
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
    poll(&mut cancelled)?;
    Ok(GasDropletTransportCandidate {
        droplets,
        parcels: next_parcels,
        catalog,
        films,
        report: GasDropletTransportReport {
            duration_s,
            gas_region,
            transport,
            coupled,
            mass_residual_kg,
            absolute_mass_defect_kg,
            energy_residual_j,
            absolute_energy_defect_j,
        },
    })
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
