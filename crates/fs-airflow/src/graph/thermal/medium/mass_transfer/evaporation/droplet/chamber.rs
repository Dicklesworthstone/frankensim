//! Original sourced spherical films sharing one finite, well-mixed atmosphere.
//!
//! Each physical droplet explicitly selects stagnant diffusion or the admitted
//! forced-sphere source adapter. All molecular and momentum sources are sampled
//! against the SAME initial gas; one original material chamber then owns their
//! shared vapor inventory, species update and gas energy. Optional reciprocal
//! sensible heat also uses one shared gas caloric solve, not successive scalar
//! updates whose result depends on the input order.
//!
//! Source claims, actual liquid calorics, finite interface kinetics and the
//! forced adapter's mean-film temperature and dilute guards remain original.
//! Initial geometry and transport freeze for this first-order interval. Every
//! final size, original gas state and authored correlation/source domain is
//! admitted before ANY candidate is returned. Endpoint samples are diagnostics
//! for the next interval, not another phase or heat transfer.
//!
//! The gas volume stays externally prescribed. This well-mixed shared humidity
//! and heat feedback does not resolve droplet wakes, shielding, spatial vapor
//! plumes, trajectories, breakup, coalescence, nucleation or complete dry-out.
//! Each row represents one actual droplet; statistical multiplicity and hidden
//! heat/mass analogies are not supplied. Original source receipts and the
//! explicitly bounded roster retain their ordinary allocation behavior.

/// Adaptive physical time with refreshed original spherical film sources.
pub mod adaptive;

/// Original gas-network transport before one joint sourced droplet interval.
pub mod transport;

use super::forced::{
    DiluteLaw, ForcedSphericalDropletFilmSample, SourceForcedSphericalDropletFilm,
};
use super::{
    ClosedGasParcel, EvaporationError, EvaporationLaw, GasFilmError, GasFilmHeatSolve,
    GasLiquidBinding, SourceSphericalDropletFilm, SphericalDropletFilmError,
    SphericalDropletFilmSample, SphericalDropletGeometry, SphericalDropletState,
    ThermochemicalGasModel, poll,
};
use fs_material::evaporation::gas::chamber::{
    GasLiquidChamberSolve, heat::GasLiquidChamberHeatSolve,
};
use fs_material::evaporation::gas::droplet::chamber::{
    GasDropletChamberCandidate, GasDropletHeatSurface, GasDropletSurface,
    advance_droplets_with_gas as advance_material,
    advance_droplets_with_gas_and_heat as advance_material_with_heat,
};

/// An explicit transport choice for one droplet, retaining its original source.
/// Forced selection includes its prescribed relative speed and dilute bound;
/// no model is selected automatically as the size or gas state changes.
#[derive(Clone, Copy, Debug)]
pub enum SourceDropletFilm<'film, 'm> {
    Stagnant(&'film SourceSphericalDropletFilm<'m>),
    Forced(SourceForcedSphericalDropletFilm<'film, 'm>),
}
impl<'film, 'm> SourceDropletFilm<'film, 'm> {
    pub const fn source(&self) -> &'film SourceSphericalDropletFilm<'m> {
        match self {
            Self::Stagnant(source) => source,
            Self::Forced(source) => source.source(),
        }
    }

    fn resolve(
        self,
        law: &dyn EvaporationLaw,
        geometry: SphericalDropletGeometry,
        surface_temperature_k: f64,
        gas: &ClosedGasParcel,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceDropletFilmSample, SphericalDropletFilmError> {
        match self {
            Self::Stagnant(source) => Ok(SourceDropletFilmSample::Stagnant(source.resolve(
                geometry,
                gas,
                &mut *cancelled,
            )?)),
            Self::Forced(source) => {
                let sample =
                    source.resolve(geometry, surface_temperature_k, gas, &mut *cancelled)?;
                source.admit_surface(law, &sample)?;
                Ok(SourceDropletFilmSample::Forced(sample))
            }
        }
    }
}

/// One physical spherical inventory with its original phase and caloric laws.
/// Stable IDs identify droplets independently of input ordering.
#[derive(Clone, Copy)]
pub struct SourceGasDropletSurface<'a, 'm> {
    pub id: u64,
    pub law: &'a dyn EvaporationLaw,
    pub binding: &'a GasLiquidBinding,
    pub droplet: SphericalDropletState,
    pub temperature_from_enthalpy: &'a dyn Fn(f64) -> Result<f64, EvaporationError>,
    pub film: SourceDropletFilm<'a, 'm>,
}

/// Independently authored sensible heat, excluding carried vapor and latent heat.
#[derive(Clone, Copy)]
pub struct SourceGasDropletHeatSurface<'a, 'm> {
    pub surface: SourceGasDropletSurface<'a, 'm>,
    pub heat_transfer_coefficient_w_m2_k: f64,
    pub heat_controls: GasFilmHeatSolve,
}

/// Full original evidence for the explicitly selected transport model.
#[derive(Clone, Debug, PartialEq)]
pub enum SourceDropletFilmSample {
    Stagnant(SphericalDropletFilmSample),
    Forced(ForcedSphericalDropletFilmSample),
}
impl SourceDropletFilmSample {
    pub const fn geometry(&self) -> SphericalDropletGeometry {
        match self {
            Self::Stagnant(sample) => sample.geometry,
            Self::Forced(sample) => sample.geometry,
        }
    }
    pub const fn binary_diffusivity_m2_s(&self) -> f64 {
        match self {
            Self::Stagnant(sample) => sample.diffusion.diffusivity_m2_s,
            Self::Forced(sample) => sample.diffusion.diffusivity_m2_s,
        }
    }
    pub const fn sherwood_number(&self) -> f64 {
        match self {
            Self::Stagnant(_) => 2.0,
            Self::Forced(sample) => sample.correlation.sherwood_number(),
        }
    }
    pub const fn mass_transfer_coefficient_m_s(&self) -> f64 {
        match self {
            Self::Stagnant(sample) => sample.mass_transfer_coefficient_m_s,
            Self::Forced(sample) => sample.mass_transfer_coefficient_m_s(),
        }
    }
}

/// Used and next-interval source samples in canonical stable-ID order.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceDropletFilmStep {
    pub id: u64,
    pub films: [SourceDropletFilmSample; 2],
}

/// All physical endpoints and source receipts are accepted together.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceGasDropletChamberCandidate<'catalog> {
    original: GasDropletChamberCandidate<'catalog>,
    films: Vec<SourceDropletFilmStep>,
}
impl<'catalog> SourceGasDropletChamberCandidate<'catalog> {
    pub const fn original(&self) -> &GasDropletChamberCandidate<'catalog> {
        &self.original
    }
    pub fn films(&self) -> &[SourceDropletFilmStep] {
        &self.films
    }
    pub fn into_parts(
        self,
    ) -> (
        GasDropletChamberCandidate<'catalog>,
        Vec<SourceDropletFilmStep>,
    ) {
        (self.original, self.films)
    }
}

/// One original joint phase interval with every film sampled at the same gas.
/// All final source/correlation support is checked atomically, including for
/// droplets whose phase interface is sealed. Empty rosters retain the original
/// material chamber's empty-interval behavior.
pub fn advance_droplets_with_gas<'catalog>(
    surfaces: &[SourceGasDropletSurface<'_, '_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasLiquidChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<SourceGasDropletChamberCandidate<'catalog>, SphericalDropletFilmError> {
    advance(
        surfaces.iter().copied().map(|surface| (surface, None)),
        gas,
        catalog,
        duration_s,
        controls,
        None,
        cancelled,
    )
}

/// One simultaneous reciprocal heat solve followed by one joint phase interval.
/// Initial geometry and sourced transport freeze across both substeps, while
/// every phase evaluation keeps the original forced-film dilute guard active.
pub fn advance_droplets_with_gas_and_heat<'catalog>(
    surfaces: &[SourceGasDropletHeatSurface<'_, '_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasLiquidChamberSolve,
    heat_controls: GasLiquidChamberHeatSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<SourceGasDropletChamberCandidate<'catalog>, SphericalDropletFilmError> {
    advance(
        surfaces.iter().map(|surface| {
            (
                surface.surface,
                Some((
                    surface.heat_transfer_coefficient_w_m2_k,
                    surface.heat_controls,
                )),
            )
        }),
        gas,
        catalog,
        duration_s,
        controls,
        Some(heat_controls),
        cancelled,
    )
}

struct Prepared<'a, 'm> {
    input: SourceGasDropletSurface<'a, 'm>,
    initial: SourceDropletFilmSample,
    guarded: Option<DiluteLaw<'a, dyn EvaporationLaw + 'a>>,
    heat: Option<(f64, GasFilmHeatSolve)>,
}
impl Prepared<'_, '_> {
    fn material_surface(&self) -> GasDropletSurface<'_> {
        let law: &dyn EvaporationLaw = self.guarded.as_ref().map_or(self.input.law, |guard| guard);
        GasDropletSurface {
            id: self.input.id,
            law,
            binding: self.input.binding,
            droplet: self.input.droplet,
            temperature_from_enthalpy: self.input.temperature_from_enthalpy,
            binary_diffusivity_m2_s: self.initial.binary_diffusivity_m2_s(),
            sherwood_number: self.initial.sherwood_number(),
        }
    }
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance<'surface, 'm: 'surface, 'catalog>(
    surfaces: impl ExactSizeIterator<
        Item = (
            SourceGasDropletSurface<'surface, 'm>,
            Option<(f64, GasFilmHeatSolve)>,
        ),
    >,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasLiquidChamberSolve,
    heat_controls: Option<GasLiquidChamberHeatSolve>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SourceGasDropletChamberCandidate<'catalog>, SphericalDropletFilmError> {
    poll(&mut cancelled)?;
    if controls.maximum_surfaces == 0 || surfaces.len() > controls.maximum_surfaces {
        return Err(GasFilmError::Invalid(
            "sourced droplet count exceeds its positive chamber work ceiling",
        )
        .into());
    }
    let mut ordered: Vec<_> = surfaces.collect();
    ordered.sort_unstable_by_key(|(surface, _)| surface.id);
    if ordered.windows(2).any(|pair| pair[0].0.id == pair[1].0.id) {
        return Err(GasFilmError::Invalid("sourced droplets require unique stable IDs").into());
    }
    let mut prepared = Vec::with_capacity(ordered.len());
    for (input, heat) in ordered {
        poll(&mut cancelled)?;
        // Identity must refuse before calling the liquid caloric source.
        input.film.source().admit_binding(input.binding)?;
        let temperature_k =
            (input.temperature_from_enthalpy)(input.droplet.specific_enthalpy_j_kg)?;
        let geometry = input.droplet.geometry(input.law, temperature_k)?;
        let initial =
            input
                .film
                .resolve(input.law, geometry, temperature_k, gas, &mut cancelled)?;
        let guarded = match (&input.film, &initial) {
            (SourceDropletFilm::Forced(film), SourceDropletFilmSample::Forced(sample)) => {
                Some(DiluteLaw {
                    source: input.law,
                    maximum: film.maximum_vapor_mole_fraction(),
                    film_temperature_k: sample.film_temperature_k,
                })
            }
            _ => None,
        };
        prepared.push(Prepared {
            input,
            initial,
            guarded,
            heat,
        });
    }
    let original = if let Some(heat_controls) = heat_controls {
        let rows: Vec<_> = prepared
            .iter()
            .map(|entry| {
                let (coefficient, heat_controls) = entry
                    .heat
                    .ok_or(GasFilmError::Invalid("sourced droplet heat row is missing"))?;
                Ok(GasDropletHeatSurface {
                    surface: entry.material_surface(),
                    heat_transfer_coefficient_w_m2_k: coefficient,
                    heat_controls,
                })
            })
            .collect::<Result<_, GasFilmError>>()?;
        advance_material_with_heat(
            &rows,
            gas,
            catalog,
            duration_s,
            controls,
            heat_controls,
            &mut cancelled,
        )?
    } else {
        let rows: Vec<_> = prepared.iter().map(Prepared::material_surface).collect();
        advance_material(&rows, gas, catalog, duration_s, controls, &mut cancelled)?
    };
    if original.droplets().len() != prepared.len() {
        return Err(GasFilmError::Invalid("original droplet chamber changed its roster").into());
    }
    let mut films = Vec::with_capacity(prepared.len());
    for (entry, endpoint) in prepared.into_iter().zip(original.droplets()) {
        poll(&mut cancelled)?;
        if entry.input.id != endpoint.id {
            return Err(
                GasFilmError::Invalid("original droplet chamber changed stable IDs").into(),
            );
        }
        let last = entry.input.film.resolve(
            entry.input.law,
            endpoint.geometry[1],
            endpoint.exchange.surface_temperature_k[1],
            original.gas(),
            &mut cancelled,
        )?;
        films.push(SourceDropletFilmStep {
            id: endpoint.id,
            films: [entry.initial, last],
        });
    }
    poll(&mut cancelled)?;
    Ok(SourceGasDropletChamberCandidate { original, films })
}
