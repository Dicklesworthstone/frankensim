//! Source-backed binary gas films for the original solid/vapor exchange owner.
//!
//! Molecular diffusivity is an explicitly selected species-pair source, never
//! thermal diffusivity or an assumed Lewis number. The original fluid supplies
//! viscosity/density at the declared representative film temperature. The
//! resulting coefficient can enter [`SurfaceKinetics`] directly; the existing
//! thermal owner keeps signed mass, inventory limits, latent heat and recession.
//!
//! The laminar flat-plate law is `Sh_L = 0.664 Re_L^(1/2) Sc^(1/3)` and
//! `k_m = Sh_L D_AB/L`. Source: Seader and Henley, *Separation Process
//! Principles*, 2nd ed., Wiley, section 3.4, Eqs. (3-132), (3-137)--(3-138):
//! <https://catalogimages.wiley.com/images/db/pdf/0471464805.ch3.pdf>.
//! The source assumes a thin attached laminar boundary layer, uniform surface
//! concentration, constant properties, zero streamwise pressure gradient and
//! sufficiently weak transfer that the velocity field is undisturbed. Its
//! transition limit is `Re_L = 5e5`; this implementation conservatively admits
//! `Sc >= 0.6`, within the source's `Sc > 0.5` analogy. Positive speed is required:
//! a stagnant gas needs a diffusion geometry, not a zero-film extrapolation.
//!
//! A strip retains the physical leading edge. Its mean is the difference of
//! cumulative transfer divided by strip length, evaluated without subtraction:
//! `k_[a,b] = k_[0,b] / (1 + sqrt(a/b))`. Thus subdivision preserves total
//! transfer for one uniform wall, fluid state and driving concentration. This
//! is a frozen-property, dilute binary continuum approximation. It does not
//! model Stefan flow, turbulence, buoyancy, nucleation or a resolved plume.

/// Live binary films inside the original finite-gas, two-temperature histories.
pub mod chamber;

/// Original laminar plate/strip physics with explicitly authored transport.
pub mod laminar_plate;

/// Source-resolved liquid evaporation on the original coupled thermal timeline.
pub mod evaporation;

use fs_matdb::{MatDbError, QueryPoint};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransportState, ThermalFluidTransport};
use fs_material::state_point::{
    MaterialStatePointError, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::{
    BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial,
    BinaryDiffusionSample,
};
use fs_qty::chemistry::SpeciesId;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Length, Pressure, QuantitySpec, Temperature};
use fs_thermal::sublimation::SurfaceKinetics;

/// Explicit geometry and representative state, held over a phase-transfer step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarPlateFilmConditions {
    /// Start measured downstream from the physical flow-facing leading edge.
    /// Zero evaluates the usual whole-plate mean.
    pub start: Length,
    /// End measured from the SAME leading edge, strictly greater than `start`.
    pub end: Length,
    /// Positive relative free-stream speed parallel to the plate [m/s].
    pub speed_m_s: f64,
    /// Actual temperature of the supplied fluid transport evaluation.
    pub temperature: Temperature,
    /// Total gas pressure used to query binary diffusivity, not vapor pressure.
    pub total_pressure: Pressure,
}

/// Source/domain refusal returns no partial film or modified physical state.
#[derive(Debug)]
pub enum GasFilmError {
    Invalid(&'static str),
    Material(MaterialError),
    /// Original coating-temperature inversion failure in a live film resolver.
    Surface(fs_thermal::sublimation::SublimationError),
    Query(MatDbError),
    Source(MaterialStatePointError),
    /// Actual volatile loading leaves the live adapter's authored dilute domain.
    OutsideDilute {
        bulk_vapor_mole_fraction: f64,
        saturation_pressure_fraction: f64,
        kinetic_equilibrium_pressure_fraction: f64,
        maximum_vapor_mole_fraction: f64,
    },
    /// No extrapolation across transition or below the admitted Schmidt range.
    OutsideCorrelation {
        reynolds: f64,
        schmidt: f64,
    },
    Cancelled,
}
impl core::fmt::Display for GasFilmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source gas film: {self:?}")
    }
}
impl std::error::Error for GasFilmError {}
impl From<laminar_plate::LaminarPlateFilmError> for GasFilmError {
    fn from(error: laminar_plate::LaminarPlateFilmError) -> Self {
        match error {
            laminar_plate::LaminarPlateFilmError::Invalid(reason) => Self::Invalid(reason),
            laminar_plate::LaminarPlateFilmError::OutsideCorrelation { reynolds, schmidt } => {
                Self::OutsideCorrelation { reynolds, schmidt }
            }
        }
    }
}

/// One material-preserving film evaluation, not a new gas or deposition model.
#[derive(Debug)]
pub struct SourceLaminarPlateGasFilm<'f, 'm, F: ThermalFluidTransport + ?Sized> {
    fluid: &'f F,
    source: BinaryDiffusionMaterial<'m>,
    diffusion: BinaryDiffusionSample,
    vapor_index: usize,
    transport: FluidTransportState,
    conditions: LaminarPlateFilmConditions,
    correlation: laminar_plate::LaminarPlateFilm,
}
impl<'m, F: ThermalFluidTransport + ?Sized> SourceLaminarPlateGasFilm<'_, 'm, F> {
    /// Original fluid owner, retaining its own material evidence and choices.
    pub const fn fluid(&self) -> &F {
        self.fluid
    }
    /// Original selected pair binding, including the unmodified query context.
    pub const fn diffusion_source(&self) -> &BinaryDiffusionMaterial<'m> {
        &self.source
    }
    /// Canonical typed species pair, actual diffusivity and complete source receipt.
    pub const fn diffusion(&self) -> &BinaryDiffusionSample {
        &self.diffusion
    }
    pub fn vapor_species(&self) -> &SpeciesId {
        &self.diffusion.species[self.vapor_index]
    }
    pub fn carrier_species(&self) -> &SpeciesId {
        &self.diffusion.species[1 - self.vapor_index]
    }
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }
    pub const fn conditions(&self) -> LaminarPlateFilmConditions {
        self.conditions
    }
    /// The same admitted strip law used by the source-free physical evaluator.
    pub const fn correlation(&self) -> &laminar_plate::LaminarPlateFilm {
        &self.correlation
    }
    /// Reynolds number at the downstream endpoint, including the upstream plate.
    pub const fn reynolds(&self) -> f64 {
        self.correlation.reynolds()
    }
    /// Actual molecular Schmidt number `nu/D_AB`, independent of thermal transport.
    pub const fn schmidt(&self) -> f64 {
        self.correlation.schmidt()
    }
    /// Mean Sherwood number over `[0, end]`, not based on a cell's own length.
    pub const fn downstream_mean_sherwood(&self) -> f64 {
        self.correlation.downstream_mean_sherwood()
    }
    /// Area-mean coefficient on the requested strip [m/s].
    pub const fn mass_transfer_coefficient_m_s(&self) -> f64 {
        self.correlation.mass_transfer_coefficient_m_s()
    }

    /// Bind this admitted film to the existing signed solid/vapor kinetics.
    ///
    /// The accommodation probability is independently authored, never inferred
    /// from flow or diffusivity. This returns the original owner type, usable by
    /// `fs_thermal::sublimation::coating::advance_supported_surface_with_heat`
    /// and the shared chamber/timeline consumers. The film remains frozen over
    /// that call. Its total pressure is distinct from the owner's volatile
    /// partial pressure; callers must bind the same volatile species and keep
    /// the dilute interface concentrations and gas/film assumptions valid
    /// throughout the accepted step. The receiving owner's gas-temperature
    /// conversion from concentration to pressure conductance is unchanged.
    pub fn surface_kinetics(
        &self,
        accommodation_coefficient: f64,
    ) -> Result<SurfaceKinetics, GasFilmError> {
        self.correlation
            .surface_kinetics(accommodation_coefficient)
            .map_err(Into::into)
    }

    /// Re-evaluate the SAME selected diffusivity law at another physical state.
    /// Source support ends in refusal; another overlapping claim is never chosen.
    /// The supplied fluid remains its original independently re-evaluated owner.
    pub fn at_film_state<'g, G: ThermalFluidTransport + ?Sized>(
        &self,
        fluid: &'g G,
        conditions: LaminarPlateFilmConditions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SourceLaminarPlateGasFilm<'g, 'm, G>, GasFilmError> {
        evaluate(
            fluid,
            &self.source,
            self.vapor_species(),
            conditions,
            Some(&self.diffusion.material),
            &mut cancelled,
        )
    }
}

/// Resolve one original binary source and admit the laminar plate film.
///
/// The two distinct canonical species must include `vapor_species`; the other
/// is the explicitly bound carrier. Material names never identify a species.
/// Only canonical `temperature` and `pressure` source coordinates are supported.
/// They are replaced by the actual film values with explicit physical kinds;
/// other context axes remain unchanged. Typed T/p aliases and the legacy
/// untyped `T`/`p` names refuse rather than retain stale applicability values.
///
/// The supplied thermal and momentum views must agree exactly and use the
/// declared temperature. The generic fluid interface has no pressure or
/// chemical-identity accessor: matching total pressure, carrier/composition,
/// continuum validity and weak-transfer boundary-layer assumptions are caller
/// declarations. A raw authored transport state does not gain fluid source
/// evidence by being paired with sourced molecular diffusivity.
///
/// `start` and `end` follow the actual flow direction; reversed flow needs
/// distances from the opposite physical edge. Local property/concentration
/// variations require a separately admitted frozen-strip approximation.
pub fn evaluate_laminar_plate_gas_film<'f, 'm, F: ThermalFluidTransport + ?Sized>(
    fluid: &'f F,
    source: &BinaryDiffusionMaterial<'m>,
    vapor_species: &SpeciesId,
    conditions: LaminarPlateFilmConditions,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SourceLaminarPlateGasFilm<'f, 'm, F>, GasFilmError> {
    evaluate(
        fluid,
        source,
        vapor_species,
        conditions,
        None,
        &mut cancelled,
    )
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasFilmError> {
    if cancelled() {
        Err(GasFilmError::Cancelled)
    } else {
        Ok(())
    }
}

fn query(
    context: &QueryPoint,
    conditions: LaminarPlateFilmConditions,
) -> Result<QueryPoint, GasFilmError> {
    if context.axes().contains_key("T") || context.axes().contains_key("p") {
        return Err(GasFilmError::Invalid(
            "legacy T/p aliases are not canonical binary gas film coordinates",
        ));
    }
    for (axis, quantity) in context.axis_quantities() {
        let expected = if quantity.dims() == Temperature::DIMS || axis == "temperature" {
            Some(("temperature", typed(QuantityKind::AbsoluteTemperature)))
        } else if quantity.dims() == Pressure::DIMS || axis == "pressure" {
            Some(("pressure", typed(QuantityKind::Pressure)))
        } else {
            None
        };
        if let Some((name, expected)) = expected
            && (axis != name || *quantity != expected)
        {
            return Err(GasFilmError::Invalid(
                "binary gas film needs one canonical absolute temperature and total pressure coordinate",
            ));
        }
    }
    context
        .clone()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            conditions.temperature.value(),
        )
        .and_then(|point| {
            point.with_quantity(
                "pressure",
                typed(QuantityKind::Pressure),
                conditions.total_pressure.value(),
            )
        })
        .map_err(GasFilmError::Query)
}

fn evaluate<'f, 'm, F: ThermalFluidTransport + ?Sized>(
    fluid: &'f F,
    source: &BinaryDiffusionMaterial<'m>,
    vapor_species: &SpeciesId,
    conditions: LaminarPlateFilmConditions,
    previous: Option<&ResolvedMaterialStatePoint>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SourceLaminarPlateGasFilm<'f, 'm, F>, GasFilmError> {
    poll(cancelled)?;
    let geometry = laminar_plate::LaminarPlateGeometry {
        start_m: conditions.start.value(),
        end_m: conditions.end.value(),
        speed_m_s: conditions.speed_m_s,
    };
    geometry.validate()?;
    if !positive(conditions.temperature.value()) || !positive(conditions.total_pressure.value()) {
        return Err(GasFilmError::Invalid(
            "positive finite film temperature and total pressure required",
        ));
    }
    let mut species = source.species.clone();
    species.sort();
    if species[0] == species[1] {
        return Err(GasFilmError::Invalid(
            "binary diffusivity requires two distinct species",
        ));
    }
    let vapor_index = species
        .iter()
        .position(|species| species == vapor_species)
        .ok_or(GasFilmError::Invalid(
            "selected vapor is absent from the binary diffusivity pair",
        ))?;
    let thermal = fluid
        .thermal_transport_state()
        .map_err(GasFilmError::Material)?;
    let transport = fluid.transport_state().map_err(GasFilmError::Material)?;
    if thermal.temperature_k() != conditions.temperature.value() || thermal.transport() != transport
    {
        return Err(GasFilmError::Invalid(
            "fluid thermal/momentum views must agree at the declared film temperature",
        ));
    }
    let point = query(source.context, conditions)?;
    poll(cancelled)?;
    let material = if let Some(previous) = previous {
        previous.at_query_point(source.material, &point)
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
        .ok_or(GasFilmError::Invalid("resolved binary diffusivity missing"))?;
    if property.source_axes().iter().any(|(axis, dims)| {
        (*dims == Temperature::DIMS && axis != "temperature")
            || (*dims == Pressure::DIMS && axis != "pressure")
    }) {
        return Err(GasFilmError::Invalid(
            "binary diffusivity depends on a noncanonical temperature or pressure axis",
        ));
    }
    let diffusivity_m2_s = property.value_si();
    let correlation = laminar_plate::evaluate(
        geometry,
        transport.kinematic_viscosity_m2_s(),
        diffusivity_m2_s,
    )?;
    Ok(SourceLaminarPlateGasFilm {
        fluid,
        source: source.clone(),
        diffusion: BinaryDiffusionSample {
            species,
            diffusivity_m2_s,
            material,
        },
        vapor_index,
        transport,
        conditions,
        correlation,
    })
}

#[cfg(test)]
mod tests;
