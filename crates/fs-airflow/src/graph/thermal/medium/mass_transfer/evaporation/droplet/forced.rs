//! Explicit weak-transfer forced-sphere films using original gas transport.
//!
//! The selected direct mass-transfer correlation lives in `mass_transfer::sphere`.
//! This adapter resolves molecular D and original gas rho/mu at the arithmetic
//! mean of the actual surface and bulk temperatures, at the actual total gas
//! pressure and composition. The same mean is the representative temperature
//! in the original Stefan film; interface kinetics keep their separate actual
//! surface temperature. Molecular D and Sherwood retain separate receipts.
//!
//! Geometry, D, momentum transport, Sherwood and representative film temperature
//! freeze for one first-order interval. Initial and accepted endpoint source
//! support are checked, and every phase-film evaluation admits the authored
//! dilute bound. The original Stefan logarithm and interface kinetics remain
//! active in this limited scope; this is not a strong-blowing correction.
//! Area-averaged Sherwood and finite kinetics form a lumped series closure;
//! angular mass-transfer and interfacial-kinetic variation are not resolved.
//! Uniform spherical surface, dilute binary continuum gas, nearly constant
//! film properties, negligible buoyancy/compressibility/deformation and a
//! quasi-steady velocity field remain caller assumptions. No velocity, drag,
//! droplet trajectory or sensible heat coefficient is inferred.

use super::{
    BinaryDiffusionSample, ClosedGasParcel, EvaporationError, EvaporationLaw,
    GasDropletCandidate, GasEvaporationSolve, GasFilmError, GasFilmHeatSolve,
    GasLiquidBinding, SourceSphericalDropletFilm, SphericalDropletFilmError,
    SphericalDropletGeometry, SphericalDropletState, ThermochemicalGasModel,
    evaluate_at_temperature, poll,
};
use crate::graph::thermal::medium::mass_transfer::sphere::{
    self, ForcedSphereFilm, ForcedSphereFilmError, ForcedSphereGeometry,
    ForcedSphereMassTransfer,
};
use crate::graph::thermal::medium::mass_transfer::typed;
use fs_material::evaporation::gas::droplet::{
    advance_droplet_with_gas_and_heat_at_sherwood, advance_droplet_with_gas_at_sherwood,
};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_qty::semantic::QuantityKind;
use fs_thermal::evaporation::{
    EvaporationProperties, GasFilm, InterfaceKinetics, SurfaceFlux, properties,
};

/// Conservative implementation ceiling for an explicitly authored dilute bound.
/// It is an admission policy, not a universal correlation accuracy guarantee.
pub const MAX_DILUTE_SPHERE_VAPOR_MOLE_FRACTION: f64 = 0.05;

/// Original molecular and momentum transport at one actual gas/size endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ForcedSphericalDropletFilmSample {
    pub geometry: SphericalDropletGeometry,
    pub diffusion: BinaryDiffusionSample,
    pub surface_temperature_k: f64,
    pub bulk_temperature_k: f64,
    pub film_temperature_k: f64,
    pub total_pressure_pa: f64,
    pub vapor_mole_fraction: f64,
    /// Original gas material resolved at mean temperature and actual bulk p.
    /// None reuses the caller's exact gas state when mean and bulk T coincide.
    /// This is a material evaluation, not a heated parcel or energy increment.
    pub film_state: Option<ResolvedThermochemicalGasState>,
    pub transport: FluidTransportState,
    pub correlation: ForcedSphereFilm,
}
impl ForcedSphericalDropletFilmSample {
    pub const fn mass_transfer_coefficient_m_s(&self) -> f64 {
        self.correlation.mass_transfer_coefficient_m_s()
    }
}

/// Original accepted droplet/gas with initial and next-interval film receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceForcedGasDropletCandidate<'catalog> {
    original: GasDropletCandidate<'catalog>,
    films: [ForcedSphericalDropletFilmSample; 2],
}
impl<'catalog> SourceForcedGasDropletCandidate<'catalog> {
    pub const fn original(&self) -> &GasDropletCandidate<'catalog> {
        &self.original
    }
    /// Index 0 is used for this interval. Index 1 is a refreshed diagnostic,
    /// not a second transfer or a retroactive change to the accepted step.
    pub const fn films(&self) -> &[ForcedSphericalDropletFilmSample; 2] {
        &self.films
    }
    pub fn into_parts(
        self,
    ) -> (
        GasDropletCandidate<'catalog>,
        [ForcedSphericalDropletFilmSample; 2],
    ) {
        (self.original, self.films)
    }
}

/// Explicit forced-flow selection borrowing the same pinned molecular source.
///
/// Relative speed is prescribed. The caller's positive dilute bound applies to
/// bulk volatile fraction and equilibrium surface pressure/total pressure.
/// Finite kinetics may reduce actual interface pressure; this admission remains
/// conservative. Material receipts retain their ordinary allocation behavior.
#[derive(Clone, Copy, Debug)]
pub struct SourceForcedSphericalDropletFilm<'film, 'm> {
    source: &'film SourceSphericalDropletFilm<'m>,
    model: ForcedSphereMassTransfer,
    relative_speed_m_s: f64,
    maximum_vapor_mole_fraction: f64,
}
impl<'film, 'm> SourceForcedSphericalDropletFilm<'film, 'm> {
    pub fn try_new(
        source: &'film SourceSphericalDropletFilm<'m>,
        model: ForcedSphereMassTransfer,
        relative_speed_m_s: f64,
        maximum_vapor_mole_fraction: f64,
    ) -> Result<Self, SphericalDropletFilmError> {
        if !relative_speed_m_s.is_finite() || relative_speed_m_s <= 0.0 {
            return Err(GasFilmError::Invalid(
                "forced spherical film requires positive finite prescribed relative speed",
            )
            .into());
        }
        if !maximum_vapor_mole_fraction.is_finite()
            || maximum_vapor_mole_fraction <= 0.0
            || maximum_vapor_mole_fraction > MAX_DILUTE_SPHERE_VAPOR_MOLE_FRACTION
        {
            return Err(GasFilmError::Invalid(
                "forced spherical film needs an explicit dilute fraction in (0,0.05]",
            )
            .into());
        }
        Ok(Self {
            source,
            model,
            relative_speed_m_s,
            maximum_vapor_mole_fraction,
        })
    }

    pub const fn source(&self) -> &'film SourceSphericalDropletFilm<'m> {
        self.source
    }
    pub const fn model(&self) -> ForcedSphereMassTransfer {
        self.model
    }
    pub const fn relative_speed_m_s(&self) -> f64 {
        self.relative_speed_m_s
    }
    pub const fn maximum_vapor_mole_fraction(&self) -> f64 {
        self.maximum_vapor_mole_fraction
    }

    /// Original rho/mu and the SAME molecular claim at mean film T and bulk p.
    /// This geometric sample admits bulk loading; advancing a liquid also
    /// admits its equilibrium surface pressure at each original phase evaluation.
    pub fn resolve(
        &self,
        geometry: SphericalDropletGeometry,
        surface_temperature_k: f64,
        gas: &ClosedGasParcel,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedSphericalDropletFilmSample, SphericalDropletFilmError> {
        poll(&mut cancelled)?;
        if !surface_temperature_k.is_finite() || surface_temperature_k <= 0.0 {
            return Err(GasFilmError::Invalid(
                "forced sphere requires a positive actual surface temperature",
            )
            .into());
        }
        let bulk = gas.state().state();
        let film_temperature_k = 0.5 * surface_temperature_k + 0.5 * bulk.temperature;
        let molecular = evaluate_at_temperature(
            self.source.diffusion_source(),
            self.source.vapor_species(),
            gas,
            Some(self.source.selected_diffusivity()),
            film_temperature_k,
            &mut cancelled,
        )?;
        self.admit_dilute(molecular.vapor_mole_fraction, 0.0)?;
        poll(&mut cancelled)?;
        let film_state = if film_temperature_k == bulk.temperature {
            None
        } else {
            let point = gas
                .point()
                .clone()
                .with_quantity(
                    "temperature",
                    typed(QuantityKind::AbsoluteTemperature),
                    film_temperature_k,
                )
                .and_then(|point| {
                    point.with_quantity("pressure", typed(QuantityKind::Pressure), bulk.pressure)
                })
                .map_err(GasFilmError::Query)?;
            Some(gas.model().evaluate(&point).map_err(GasFilmError::Material)?)
        };
        poll(&mut cancelled)?;
        let transport = film_state
            .as_ref()
            .unwrap_or_else(|| gas.state())
            .transport_state()
            .map_err(GasFilmError::Material)?;
        let correlation = sphere::evaluate(
            self.model,
            ForcedSphereGeometry {
                diameter_m: geometry.diameter_m(),
                relative_speed_m_s: self.relative_speed_m_s,
            },
            transport.kinematic_viscosity_m2_s(),
            molecular.diffusion.diffusivity_m2_s,
        )
        .map_err(correlation_error)?;
        poll(&mut cancelled)?;
        Ok(ForcedSphericalDropletFilmSample {
            geometry,
            diffusion: molecular.diffusion,
            surface_temperature_k,
            bulk_temperature_k: bulk.temperature,
            film_temperature_k,
            total_pressure_pa: bulk.pressure,
            vapor_mole_fraction: molecular.vapor_mole_fraction,
            film_state,
            transport,
            correlation,
        })
    }

    /// Original finite-gas phase interval with explicit dilute forced transport.
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
    ) -> Result<SourceForcedGasDropletCandidate<'catalog>, SphericalDropletFilmError> {
        let initial = self.initial(law, binding, droplet, gas, &temperature, &mut cancelled)?;
        let guarded = DiluteLaw {
            source: law,
            maximum: self.maximum_vapor_mole_fraction,
            film_temperature_k: initial.film_temperature_k,
        };
        let original = advance_droplet_with_gas_at_sherwood(
            &guarded,
            binding,
            droplet,
            gas,
            catalog,
            temperature,
            initial.diffusion.diffusivity_m2_s,
            initial.correlation.sherwood_number(),
            duration_s,
            controls,
            &mut cancelled,
        )?;
        self.finish(law, original, initial, &mut cancelled)
    }

    /// Original reciprocal heat-first/phase split with independently authored h.
    /// The initial forced film is frozen across both substeps. No heat/mass
    /// analogy or additional sensible/vapor energy account is introduced.
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
    ) -> Result<SourceForcedGasDropletCandidate<'catalog>, SphericalDropletFilmError> {
        let initial = self.initial(law, binding, droplet, gas, &temperature, &mut cancelled)?;
        let guarded = DiluteLaw {
            source: law,
            maximum: self.maximum_vapor_mole_fraction,
            film_temperature_k: initial.film_temperature_k,
        };
        let original = advance_droplet_with_gas_and_heat_at_sherwood(
            &guarded,
            binding,
            droplet,
            gas,
            catalog,
            temperature,
            initial.diffusion.diffusivity_m2_s,
            initial.correlation.sherwood_number(),
            heat_transfer_coefficient_w_m2_k,
            duration_s,
            controls,
            heat_controls,
            &mut cancelled,
        )?;
        self.finish(law, original, initial, &mut cancelled)
    }

    fn initial<L: EvaporationLaw + ?Sized>(
        &self,
        law: &L,
        binding: &GasLiquidBinding,
        droplet: SphericalDropletState,
        gas: &ClosedGasParcel,
        temperature: &impl Fn(f64) -> Result<f64, EvaporationError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ForcedSphericalDropletFilmSample, SphericalDropletFilmError> {
        poll(cancelled)?;
        self.source.admit_binding(binding)?;
        let temperature_k = temperature(droplet.specific_enthalpy_j_kg)?;
        let geometry = droplet.geometry(law, temperature_k)?;
        let initial = self.resolve(geometry, temperature_k, gas, &mut *cancelled)?;
        self.admit_surface(law, &initial)?;
        poll(cancelled)?;
        Ok(initial)
    }

    fn finish<'catalog, L: EvaporationLaw + ?Sized>(
        &self,
        law: &L,
        original: GasDropletCandidate<'catalog>,
        initial: ForcedSphericalDropletFilmSample,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceForcedGasDropletCandidate<'catalog>, SphericalDropletFilmError> {
        poll(cancelled)?;
        let last = self.resolve(
            original.geometry()[1],
            original.step().surface_temperature_k[1],
            original.gas(),
            &mut *cancelled,
        )?;
        self.admit_surface(law, &last)?;
        poll(cancelled)?;
        Ok(SourceForcedGasDropletCandidate {
            original,
            films: [initial, last],
        })
    }

    fn admit_surface<L: EvaporationLaw + ?Sized>(
        &self,
        law: &L,
        film: &ForcedSphericalDropletFilmSample,
    ) -> Result<(), SphericalDropletFilmError> {
        let surface = properties(law, film.surface_temperature_k)?.saturation_pressure_pa
            / film.total_pressure_pa;
        self.admit_dilute(film.vapor_mole_fraction, surface)
    }

    fn admit_dilute(&self, bulk: f64, surface: f64) -> Result<(), SphericalDropletFilmError> {
        if !bulk.is_finite()
            || !surface.is_finite()
            || bulk > self.maximum_vapor_mole_fraction
            || surface > self.maximum_vapor_mole_fraction
        {
            return Err(GasFilmError::OutsideDilute {
                bulk_vapor_mole_fraction: bulk,
                saturation_pressure_fraction: surface,
                kinetic_equilibrium_pressure_fraction: surface,
                maximum_vapor_mole_fraction: self.maximum_vapor_mole_fraction,
            }
            .into());
        }
        Ok(())
    }
}

fn correlation_error(error: ForcedSphereFilmError) -> SphericalDropletFilmError {
    match error {
        ForcedSphereFilmError::Invalid(reason) => GasFilmError::Invalid(reason),
        ForcedSphereFilmError::OutsideCorrelation { reynolds, schmidt } => {
            GasFilmError::OutsideCorrelation { reynolds, schmidt }
        }
    }
    .into()
}

/// Admit each phase-film state and supply the frozen mean transport temperature.
/// All constitutive values, actual interface temperature, kinetics and hooks
/// remain original; gas and carried-vapor calorics are untouched.
struct DiluteLaw<'a, L: EvaporationLaw + ?Sized> {
    source: &'a L,
    maximum: f64,
    film_temperature_k: f64,
}
impl<L: EvaporationLaw + ?Sized> EvaporationLaw for DiluteLaw<'_, L> {
    fn temperature_range_k(&self) -> [f64; 2] {
        self.source.temperature_range_k()
    }
    fn properties(&self, temperature_k: f64) -> Result<EvaporationProperties, EvaporationError> {
        self.source.properties(temperature_k)
    }
    fn interface_kinetics(&self) -> InterfaceKinetics {
        self.source.interface_kinetics()
    }
    fn evaluate_film(
        &self,
        admitted_properties: EvaporationProperties,
        surface_temperature_k: f64,
        mut film: GasFilm,
        kinetics: InterfaceKinetics,
    ) -> Result<SurfaceFlux, EvaporationError> {
        let bulk = film.vapor_partial_pressure_pa / film.total_pressure_pa;
        let surface = admitted_properties.saturation_pressure_pa / film.total_pressure_pa;
        if !bulk.is_finite() || !surface.is_finite()
            || bulk > self.maximum || surface > self.maximum
        {
            return Err(EvaporationError::Invalid(
                "forced sphere phase evaluation leaves the authored dilute vapor domain",
            ));
        }
        film.temperature_k = self.film_temperature_k;
        self.source.evaluate_film(admitted_properties, surface_temperature_k, film, kinetics)
    }
}
