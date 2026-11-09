//! Original molecular transport inside finite-gas liquid evaporation solves.
//!
//! The existing laminar plate/strip owner supplies `k_m` from the selected
//! binary diffusivity, actual gas viscosity/density and authored relative
//! speed. Every nonlinear and adaptive trial refreshes those properties at
//! its own gas composition, total pressure and explicit film temperature.
//! The original thermal/material owners still solve phase flux, liquid heat,
//! finite vapor inventory and gas calorics together; this adds no heat source.
//!
//! The boundary-layer correlation requires a planar attached laminar flow,
//! dilute vapor and weak normal mass transfer. The existing five-percent
//! admission cap is conservative policy, not an error certificate. The
//! underlying phase owner retains its Stefan correction; this does not extend
//! the plate correlation to strong blowing, droplets, turbulence or plumes.

use std::cell::RefCell;
use std::collections::BTreeMap;

use super::chamber::{
    LaminarPlateFilmBinding, MAX_DILUTE_VAPOR_MOLE_FRACTION, binary_vapor_fraction,
    evaluate_original_film,
};
use super::{GasFilmError, ResolvedMaterialStatePoint};
use fs_material::evaporation::gas::implicit::{
    ImplicitGasCandidate, ImplicitGasControls, ImplicitGasError, ImplicitGasSurface,
    adaptive::{
        AdaptiveImplicitGasCandidate, AdaptiveImplicitGasControls, AdaptiveImplicitGasError,
    },
    solve_deposition_with_gas_and_films,
    timeline::advance_deposition_with_gas_timeline_and_films,
};
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_thermal::enthalpy::backward_euler::BoundaryExchangeStepCoefficients;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use fs_thermal::evaporation::{EvaporationError, GasFilm, LiquidSurfaceState};

/// Source refusals retain their original cause rather than becoming a failed
/// nonlinear iteration or a substitute prescribed transfer coefficient.
#[derive(Debug)]
pub enum EvaporationFilmError {
    Film(GasFilmError),
    Liquid(EvaporationError),
    Thermal(DepositionError),
    Coupled(ImplicitGasError),
    Adaptive(AdaptiveImplicitGasError),
    OutsideDilute {
        bulk_vapor_mole_fraction: f64,
        saturation_pressure_fraction: f64,
        maximum_vapor_mole_fraction: f64,
    },
}
impl core::fmt::Display for EvaporationFilmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-resolved liquid gas film: {self:?}")
    }
}
impl std::error::Error for EvaporationFilmError {}
impl From<GasFilmError> for EvaporationFilmError {
    fn from(error: GasFilmError) -> Self {
        Self::Film(error)
    }
}

struct SelectedFilm<'source> {
    id: CellId,
    surface_index: usize,
    binding: LaminarPlateFilmBinding<'source>,
    selected: ResolvedMaterialStatePoint,
}

/// Immutable source choices and original liquid declarations for one thermal
/// graph. Bindings may be sparse; unbound interfaces retain their authored
/// coefficient and the original gas owner's representative-temperature rule.
pub struct SourceLaminarPlateEvaporationFilms<'liquid, 'source> {
    surfaces: Vec<ImplicitGasSurface<'liquid>>,
    films: Vec<SelectedFilm<'source>>,
    maximum_surfaces: usize,
}

impl<'liquid, 'source> SourceLaminarPlateEvaporationFilms<'liquid, 'source> {
    /// Pin each selected diffusivity claim once at the actual initial state.
    ///
    /// Film-bound cells must already be active. Future unbound liquid births
    /// can remain in `surfaces` for the original event timeline. All source
    /// contexts must use only canonical temperature/pressure coordinates;
    /// evolving composition cannot leave another state axis silently frozen.
    /// A carrier-only gas needs no numerical seed of the volatile species.
    /// Geometry is a caller-declared planar strip with the physical leading
    /// edge; exposed area never invents a new streamwise length.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bindings: &BTreeMap<CellId, LaminarPlateFilmBinding<'source>>,
        surfaces: &[ImplicitGasSurface<'liquid>],
        thermal: &DepositionSimulation<'_>,
        gas: &ClosedGasParcel,
        maximum_surfaces: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, EvaporationFilmError> {
        poll(&mut cancelled)?;
        if surfaces.len() > maximum_surfaces || bindings.len() > surfaces.len() {
            return Err(invalid(
                "liquid film roster exceeds its declared surface bound",
            ));
        }
        let mut ordered = BTreeMap::new();
        for surface in surfaces {
            poll(&mut cancelled)?;
            if ordered.insert(surface.cell_id, *surface).is_some() {
                return Err(invalid("duplicate liquid film surface identity"));
            }
        }
        let surfaces: Vec<_> = ordered.into_values().collect();
        let mut films = Vec::with_capacity(bindings.len());
        for (&id, binding) in bindings {
            poll(&mut cancelled)?;
            let surface_index = surfaces
                .binary_search_by_key(&id, |surface| surface.cell_id)
                .map_err(|_| invalid("liquid film binding names an unknown surface"))?;
            if binding
                .source
                .context
                .axes()
                .keys()
                .any(|axis| axis != "temperature" && axis != "pressure")
            {
                return Err(invalid(
                    "live binary diffusivity requires only canonical temperature/pressure context axes",
                ));
            }
            let (_, _, selected) = evaluate(
                binding,
                &surfaces[surface_index],
                thermal,
                gas,
                None,
                &mut cancelled,
            )?;
            films.push(SelectedFilm {
                id,
                surface_index,
                binding: binding.clone(),
                selected,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            surfaces,
            films,
            maximum_surfaces,
        })
    }

    /// Original surface/source declarations in canonical native cell order.
    pub fn surfaces(&self) -> &[ImplicitGasSurface<'liquid>] {
        &self.surfaces
    }

    /// The pinned diffusivity source and initial applicability receipt.
    pub fn selected_diffusivity(&self, id: CellId) -> Option<&ResolvedMaterialStatePoint> {
        self.films
            .binary_search_by_key(&id, |film| film.id)
            .ok()
            .map(|index| &self.films[index].selected)
    }

    /// Retained original molecular source, geometry and property convention.
    pub fn binding(&self, id: CellId) -> Option<&LaminarPlateFilmBinding<'source>> {
        self.films
            .binary_search_by_key(&id, |film| film.id)
            .ok()
            .map(|index| &self.films[index].binding)
    }

    /// Refresh each bound film from the SAME nonlinear thermal/gas endpoint.
    /// Only representative temperature and `k_m` change. Total and volatile
    /// pressures remain the original finite gas owner's exact values.
    ///
    /// The film map is speculative scratch. Any refusal discards it; neither
    /// original physical inventory nor prepared source selection is mutated.
    /// The coupled owner additionally verifies the complete active roster and
    /// original gas pressures after this callback returns.
    pub fn resolve(
        &self,
        thermal: &DepositionSimulation<'_>,
        gas: &ClosedGasParcel,
        films: &mut BTreeMap<CellId, GasFilm>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), EvaporationFilmError> {
        poll(&mut cancelled)?;
        if films.len() > self.maximum_surfaces {
            return Err(invalid(
                "liquid film roster exceeds its declared surface bound",
            ));
        }
        for bound in &self.films {
            poll(&mut cancelled)?;
            let film = films
                .get_mut(&bound.id)
                .ok_or_else(|| invalid("prepared liquid film is absent from the current roster"))?;
            let (temperature, coefficient, _) = evaluate(
                &bound.binding,
                &self.surfaces[bound.surface_index],
                thermal,
                gas,
                Some(&bound.selected),
                &mut cancelled,
            )?;
            film.temperature_k = temperature;
            film.mass_transfer_coefficient_m_s = coefficient;
        }
        poll(&mut cancelled)
    }

    /// Solve one joint liquid/thermal/gas endpoint with live molecular films.
    /// Original heat, latent energy, upstream vapor enthalpy and mass accounts
    /// remain with the finite-gas owner. `resolve_heat` supplies independent
    /// heat sources; it must not add phase heat or the surface's sensible film
    /// coefficient again. No speculative callback may advance accepted history.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn solve<'m, 'catalog, T>(
        &self,
        thermal: &DepositionSimulation<'m>,
        gas: &ClosedGasParcel,
        catalog: &'catalog ThermochemicalGasModel,
        target_s: f64,
        controls: ImplicitGasControls,
        resolve_heat: impl FnMut(
            &DepositionSimulation<'m>,
            &BTreeMap<CellId, LiquidSurfaceState>,
            &ClosedGasParcel,
        )
            -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ImplicitGasCandidate<'m, 'catalog, T>, EvaporationFilmError> {
        let cancellation = RefCell::new(&mut cancelled);
        let failure = RefCell::new(None);
        let result = solve_deposition_with_gas_and_films(
            thermal,
            &self.surfaces,
            gas,
            catalog,
            target_s,
            controls,
            resolve_heat,
            |thermal, gas, films| {
                self.resolve(thermal, gas, films, || (cancellation.borrow_mut())())
                    .map_err(|error| {
                        *failure.borrow_mut() = Some(error);
                        ImplicitGasError::Invalid("original liquid gas-film source refused")
                    })
            },
            || checkpoint(&mut || (cancellation.borrow_mut())()),
        );
        match failure.into_inner() {
            Some(error) => Err(error),
            None => result.map_err(EvaporationFilmError::Coupled),
        }
    }

    /// Advance physical time and native events with the same source feedback.
    ///
    /// Coarse, fine and rejected trials all re-evaluate the pinned molecular
    /// sources at their own states; only accepted fine endpoints transfer mass
    /// and heat. The original timeline controls physical local error and clips
    /// at authored events. Bound interfaces must stay active; unbound births
    /// retain native event semantics. The entire request is transactional.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance<'m, 'catalog, T>(
        &self,
        thermal: &DepositionSimulation<'m>,
        gas: &ClosedGasParcel,
        catalog: &'catalog ThermochemicalGasModel,
        target_s: f64,
        controls: AdaptiveImplicitGasControls,
        resolve_heat: impl FnMut(
            &DepositionSimulation<'m>,
            &BTreeMap<CellId, LiquidSurfaceState>,
            &ClosedGasParcel,
        )
            -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveImplicitGasCandidate<'m, 'catalog, T>, EvaporationFilmError> {
        let cancellation = RefCell::new(&mut cancelled);
        let failure = RefCell::new(None);
        let result = advance_deposition_with_gas_timeline_and_films(
            thermal,
            &self.surfaces,
            gas,
            catalog,
            target_s,
            controls,
            resolve_heat,
            |thermal, gas, films| {
                self.resolve(thermal, gas, films, || (cancellation.borrow_mut())())
                    .map_err(|error| {
                        *failure.borrow_mut() = Some(error);
                        ImplicitGasError::Invalid("original liquid gas-film source refused")
                    })
            },
            || checkpoint(&mut || (cancellation.borrow_mut())()),
        );
        match failure.into_inner() {
            Some(error) => Err(error),
            None => result.map_err(EvaporationFilmError::Adaptive),
        }
    }
}

fn invalid(message: &'static str) -> EvaporationFilmError {
    GasFilmError::Invalid(message).into()
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), EvaporationFilmError> {
    if cancelled() {
        Err(GasFilmError::Cancelled.into())
    } else {
        Ok(())
    }
}

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), DepositionError> {
    if cancelled() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}

fn evaluate(
    binding: &LaminarPlateFilmBinding<'_>,
    surface: &ImplicitGasSurface<'_>,
    thermal: &DepositionSimulation<'_>,
    gas: &ClosedGasParcel,
    previous: Option<&ResolvedMaterialStatePoint>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64, ResolvedMaterialStatePoint), EvaporationFilmError> {
    poll(cancelled)?;
    let ceiling = binding.maximum_vapor_mole_fraction;
    if !ceiling.is_finite() || ceiling <= 0.0 || ceiling > MAX_DILUTE_VAPOR_MOLE_FRACTION {
        return Err(invalid(
            "live dilute film needs an authored vapor-fraction ceiling in (0,0.05]",
        ));
    }
    let vapor = &surface.binding.species;
    let vapor_fraction = binary_vapor_fraction(&binding.source, vapor, gas, cancelled)?;
    let temperature = thermal
        .temperature_k(surface.cell_id)
        .map_err(EvaporationFilmError::Thermal)?;
    poll(cancelled)?;
    let phase = fs_thermal::evaporation::properties(surface.law, temperature)
        .map_err(EvaporationFilmError::Liquid)?;
    let saturation_fraction = phase.saturation_pressure_pa / gas.state().state().pressure;
    if [vapor_fraction, saturation_fraction]
        .into_iter()
        .any(|fraction| !fraction.is_finite() || fraction < 0.0 || fraction > ceiling)
    {
        return Err(EvaporationFilmError::OutsideDilute {
            bulk_vapor_mole_fraction: vapor_fraction,
            saturation_pressure_fraction: saturation_fraction,
            maximum_vapor_mole_fraction: ceiling,
        });
    }
    evaluate_original_film(binding, vapor, temperature, gas, previous, cancelled)
        .map_err(Into::into)
}
