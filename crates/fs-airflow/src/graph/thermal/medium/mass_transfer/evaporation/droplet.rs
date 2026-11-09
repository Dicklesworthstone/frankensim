//! Source-backed stagnant spherical films on the original finite-gas owner.
//!
//! The molecular source is one explicitly identified species pair and one
//! immutable selected diffusivity claim. Each call samples that SAME claim at
//! the actual bulk gas temperature and TOTAL pressure, then the original
//! spherical geometry supplies `k_m = D_AB/r` (diameter-based `Sh = 2`). This
//! exact stagnant-sphere reduction needs no viscosity, thermal diffusivity,
//! Lewis-number assumption or empirical forced-convection coefficient.
//!
//! The material droplet owner retains the original liquid enthalpy/density,
//! finite gas, nonlinear Stefan/interface kinetics and all conserved accounts.
//! Its initial area and D/r remain frozen through the whole first-order
//! interval, including optional heat-then-phase splitting. The accepted final
//! gas and size are checked against the same molecular source. That final
//! sample is a diagnostic for the next interval, not a second transfer.
//!
//! Geometry evaluation uses the original allocation-free scalar owner. Source
//! selection and samples retain the existing material resolver's receipts and
//! allocation behavior; no additional film map or patch roster is created.
//! The isolated-sphere, uniform-liquid, binary continuum assumptions and
//! exclusions in `fs_thermal::evaporation::droplet` remain authoritative.

/// Explicit dilute forced convection with the same original molecular source.
pub mod forced;

use super::super::chamber::binary_vapor_fraction;
use super::super::{GasFilmError, typed};
use fs_matdb::QueryPoint;
use fs_material::evaporation::gas::droplet::{
    GasDropletCandidate, SphericalDropletGeometry, SphericalDropletState, advance_droplet_with_gas,
    advance_droplet_with_gas_and_heat,
};
use fs_material::evaporation::gas::{
    GasEvaporationError, GasEvaporationSolve, GasFilmHeatSolve, GasLiquidBinding,
};
use fs_material::state_point::{
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::{
    BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial,
    BinaryDiffusionSample,
};
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_qty::chemistry::SpeciesId;
use fs_qty::semantic::QuantityKind;
use fs_qty::{Pressure, Temperature};
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw};

/// Original source and phase refusals retain their precise owning error.
#[derive(Debug)]
pub enum SphericalDropletFilmError {
    Film(GasFilmError),
    Liquid(EvaporationError),
    Phase(GasEvaporationError),
}
impl core::fmt::Display for SphericalDropletFilmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-resolved spherical droplet film: {self:?}")
    }
}
impl std::error::Error for SphericalDropletFilmError {}
impl From<GasFilmError> for SphericalDropletFilmError {
    fn from(error: GasFilmError) -> Self {
        Self::Film(error)
    }
}
impl From<EvaporationError> for SphericalDropletFilmError {
    fn from(error: EvaporationError) -> Self {
        Self::Liquid(error)
    }
}
impl From<GasEvaporationError> for SphericalDropletFilmError {
    fn from(error: GasEvaporationError) -> Self {
        Self::Phase(error)
    }
}

/// Original source evidence at one actual gas state and accepted geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct SphericalDropletFilmSample {
    pub geometry: SphericalDropletGeometry,
    pub diffusion: BinaryDiffusionSample,
    pub temperature_k: f64,
    pub total_pressure_pa: f64,
    pub vapor_mole_fraction: f64,
    pub mass_transfer_coefficient_m_s: f64,
}

/// Original physical candidate plus its initial and final molecular samples.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceGasDropletCandidate<'catalog> {
    original: GasDropletCandidate<'catalog>,
    films: [SphericalDropletFilmSample; 2],
}
impl<'catalog> SourceGasDropletCandidate<'catalog> {
    pub const fn original(&self) -> &GasDropletCandidate<'catalog> {
        &self.original
    }
    /// Index 0 supplies this interval's D/r. Index 1 is the actual accepted
    /// endpoint's refreshed coefficient, ready for the next physical interval.
    pub const fn films(&self) -> &[SphericalDropletFilmSample; 2] {
        &self.films
    }
    pub fn into_parts(
        self,
    ) -> (
        GasDropletCandidate<'catalog>,
        [SphericalDropletFilmSample; 2],
    ) {
        (self.original, self.films)
    }
}

/// One pinned molecular source, independent of any particular droplet size.
///
/// Only canonical temperature/pressure source coordinates are admitted. A
/// changing composition or other state axis cannot remain silently frozen.
/// Current active gas constituents must belong to the authored pair and retain
/// positive carrier inventory. An initially absent vapor needs no artificial
/// seed. The full Stefan sphere has no inherited dilute flat-plate cutoff.
#[derive(Debug)]
pub struct SourceSphericalDropletFilm<'m> {
    source: BinaryDiffusionMaterial<'m>,
    vapor_index: usize,
    selected: ResolvedMaterialStatePoint,
}
impl<'m> SourceSphericalDropletFilm<'m> {
    /// Pin one admissible original diffusivity claim at the actual initial gas.
    pub fn new(
        source: &BinaryDiffusionMaterial<'m>,
        vapor_species: &SpeciesId,
        gas: &ClosedGasParcel,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SphericalDropletFilmError> {
        let sample = evaluate(source, vapor_species, gas, None, &mut cancelled)?;
        let vapor_index = source
            .species
            .iter()
            .position(|species| species == vapor_species)
            .ok_or(GasFilmError::Invalid(
                "spherical film source omits its volatile species",
            ))?;
        poll(&mut cancelled)?;
        Ok(Self {
            source: source.clone(),
            vapor_index,
            selected: sample.diffusion.material,
        })
    }

    pub const fn diffusion_source(&self) -> &BinaryDiffusionMaterial<'m> {
        &self.source
    }
    pub const fn selected_diffusivity(&self) -> &ResolvedMaterialStatePoint {
        &self.selected
    }
    pub fn vapor_species(&self) -> &SpeciesId {
        &self.source.species[self.vapor_index]
    }
    pub fn carrier_species(&self) -> &SpeciesId {
        &self.source.species[1 - self.vapor_index]
    }

    /// Evaluate D/r at this actual gas state without changing source selection.
    pub fn resolve(
        &self,
        geometry: SphericalDropletGeometry,
        gas: &ClosedGasParcel,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SphericalDropletFilmSample, SphericalDropletFilmError> {
        let sample = evaluate(
            &self.source,
            self.vapor_species(),
            gas,
            Some(&self.selected),
            &mut cancelled,
        )?;
        let result = sample.at_geometry(geometry)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn admit_binding(&self, binding: &GasLiquidBinding) -> Result<(), SphericalDropletFilmError> {
        if &binding.species != self.vapor_species() {
            return Err(GasFilmError::Invalid(
                "spherical film cannot change its pinned volatile species",
            )
            .into());
        }
        Ok(())
    }

    /// Original finite-gas phase transfer with molecular D sampled at actual T,p.
    /// Source support, species identity, cancellation and late endpoint refusal
    /// discard the whole private candidate, leaving both caller inputs intact.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_gas<'catalog, L: EvaporationLaw + ?Sized>(
        &self,
        law: &L,
        binding: &GasLiquidBinding,
        droplet: SphericalDropletState,
        gas: &ClosedGasParcel,
        catalog: &'catalog ThermochemicalGasModel,
        temperature: impl Fn(f64) -> Result<f64, EvaporationError>,
        duration_s: f64,
        controls: GasEvaporationSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SourceGasDropletCandidate<'catalog>, SphericalDropletFilmError> {
        poll(&mut cancelled)?;
        self.admit_binding(binding)?;
        let initial = evaluate(
            &self.source,
            self.vapor_species(),
            gas,
            Some(&self.selected),
            &mut cancelled,
        )?;
        let original = advance_droplet_with_gas(
            law,
            binding,
            droplet,
            gas,
            catalog,
            temperature,
            initial.diffusion.diffusivity_m2_s,
            duration_s,
            controls,
            &mut cancelled,
        )?;
        self.finish(original, initial, &mut cancelled)
    }

    /// Delegate the original reciprocal sensible-heat then phase interval.
    /// Initial D/r is frozen across both substeps, as in the material owner.
    /// The sensible coefficient is explicit and excludes carried vapor/latent
    /// energy; this adapter does not infer it from molecular diffusivity.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_gas_and_heat<'catalog, L: EvaporationLaw + ?Sized>(
        &self,
        law: &L,
        binding: &GasLiquidBinding,
        droplet: SphericalDropletState,
        gas: &ClosedGasParcel,
        catalog: &'catalog ThermochemicalGasModel,
        temperature: impl Fn(f64) -> Result<f64, EvaporationError>,
        heat_transfer_coefficient_w_m2_k: f64,
        duration_s: f64,
        controls: GasEvaporationSolve,
        heat_controls: GasFilmHeatSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SourceGasDropletCandidate<'catalog>, SphericalDropletFilmError> {
        poll(&mut cancelled)?;
        self.admit_binding(binding)?;
        let initial = evaluate(
            &self.source,
            self.vapor_species(),
            gas,
            Some(&self.selected),
            &mut cancelled,
        )?;
        let original = advance_droplet_with_gas_and_heat(
            law,
            binding,
            droplet,
            gas,
            catalog,
            temperature,
            initial.diffusion.diffusivity_m2_s,
            heat_transfer_coefficient_w_m2_k,
            duration_s,
            controls,
            heat_controls,
            &mut cancelled,
        )?;
        self.finish(original, initial, &mut cancelled)
    }

    fn finish<'catalog>(
        &self,
        original: GasDropletCandidate<'catalog>,
        initial: DiffusionState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceGasDropletCandidate<'catalog>, SphericalDropletFilmError> {
        poll(cancelled)?;
        let first = initial.at_geometry(original.geometry()[0])?;
        let last = self.resolve(original.geometry()[1], original.gas(), &mut *cancelled)?;
        poll(cancelled)?;
        Ok(SourceGasDropletCandidate {
            original,
            films: [first, last],
        })
    }
}

struct DiffusionState {
    diffusion: BinaryDiffusionSample,
    temperature_k: f64,
    total_pressure_pa: f64,
    vapor_mole_fraction: f64,
}
impl DiffusionState {
    fn at_geometry(
        self,
        geometry: SphericalDropletGeometry,
    ) -> Result<SphericalDropletFilmSample, SphericalDropletFilmError> {
        let coefficient =
            geometry.mass_transfer_coefficient_m_s(self.diffusion.diffusivity_m2_s)?;
        Ok(SphericalDropletFilmSample {
            geometry,
            diffusion: self.diffusion,
            temperature_k: self.temperature_k,
            total_pressure_pa: self.total_pressure_pa,
            vapor_mole_fraction: self.vapor_mole_fraction,
            mass_transfer_coefficient_m_s: coefficient,
        })
    }
}

fn poll(cancelled: &mut dyn FnMut() -> bool) -> Result<(), SphericalDropletFilmError> {
    if cancelled() {
        Err(GasFilmError::Cancelled.into())
    } else {
        Ok(())
    }
}

fn query(
    source: &BinaryDiffusionMaterial<'_>,
    temperature: f64,
    pressure: f64,
) -> Result<QueryPoint, GasFilmError> {
    if source
        .context
        .axes()
        .keys()
        .any(|axis| axis != "temperature" && axis != "pressure")
    {
        return Err(GasFilmError::Invalid(
            "live spherical diffusivity requires only canonical temperature/pressure context axes",
        ));
    }
    for (axis, quantity) in source.context.axis_quantities() {
        let expected = if axis == "temperature" {
            typed(QuantityKind::AbsoluteTemperature)
        } else {
            typed(QuantityKind::Pressure)
        };
        if *quantity != expected {
            return Err(GasFilmError::Invalid(
                "spherical diffusivity needs typed absolute temperature and total pressure",
            ));
        }
    }
    source
        .context
        .clone()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .and_then(|point| point.with_quantity("pressure", typed(QuantityKind::Pressure), pressure))
        .map_err(GasFilmError::Query)
}

fn evaluate(
    source: &BinaryDiffusionMaterial<'_>,
    vapor: &SpeciesId,
    gas: &ClosedGasParcel,
    previous: Option<&ResolvedMaterialStatePoint>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<DiffusionState, SphericalDropletFilmError> {
    evaluate_at_temperature(
        source, vapor, gas, previous, gas.state().state().temperature, cancelled,
    )
}

/// Query the same admitted molecular source at an explicit transport temperature.
/// The original gas supplies composition and total pressure without heat/storage.
fn evaluate_at_temperature(
    source: &BinaryDiffusionMaterial<'_>,
    vapor: &SpeciesId,
    gas: &ClosedGasParcel,
    previous: Option<&ResolvedMaterialStatePoint>,
    temperature_k: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<DiffusionState, SphericalDropletFilmError> {
    poll(cancelled)?;
    let vapor_mole_fraction = binary_vapor_fraction(source, vapor, gas, cancelled)?;
    let bulk = gas.state().state();
    if !bulk.temperature.is_finite()
        || bulk.temperature <= 0.0
        || !bulk.pressure.is_finite()
        || bulk.pressure <= 0.0
        || vapor_mole_fraction >= 1.0
    {
        return Err(GasFilmError::Invalid(
            "spherical film requires positive finite gas temperature/pressure and a represented carrier fraction",
        ).into());
    }
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err(GasFilmError::Invalid(
            "spherical molecular source needs a positive finite transport temperature",
        ).into());
    }
    let point = query(source, temperature_k, bulk.pressure)?;
    poll(cancelled)?;
    let material = if let Some(selected) = previous {
        selected.at_query_point(source.material, &point)
    } else {
        let requirement = ScalarPropertyRequirement::try_new(
            BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .map_err(GasFilmError::Source)?;
        resolve_material_state_point(
            source.material,
            &point,
            &[requirement],
            source.selection.clone(),
        )
    }
    .map_err(GasFilmError::Source)?;
    poll(cancelled)?;
    let property = material
        .property(BINARY_DIFFUSIVITY_PROPERTY)
        .ok_or(GasFilmError::Invalid(
            "resolved spherical binary diffusivity is missing",
        ))?;
    if property.source_axes().iter().any(|(axis, dims)| {
        (*dims == Temperature::DIMS && axis != "temperature")
            || (*dims == Pressure::DIMS && axis != "pressure")
    }) {
        return Err(GasFilmError::Invalid(
            "spherical diffusivity depends on a noncanonical temperature or pressure axis",
        )
        .into());
    }
    let diffusivity_m2_s = property.value_si();
    if !diffusivity_m2_s.is_finite() || diffusivity_m2_s <= 0.0 {
        return Err(GasFilmError::Invalid("positive finite molecular diffusivity required").into());
    }
    let mut species = source.species.clone();
    species.sort();
    poll(cancelled)?;
    Ok(DiffusionState {
        diffusion: BinaryDiffusionSample {
            species,
            diffusivity_m2_s,
            material,
        },
        temperature_k,
        total_pressure_pa: bulk.pressure,
        vapor_mole_fraction,
    })
}
