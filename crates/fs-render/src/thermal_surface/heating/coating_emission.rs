//! Accepted native states determine reciprocal, nonisothermal coating emission.
//!
//! Preparation resolves the same optical laws and receiving cells as optical
//! heating. Emission does not depend on an illumination schedule. A snapshot
//! owns its resolved optics and temperatures, so later directional/spectral
//! queries do not read a changing simulation or repeat constitutive resolution.
//! Observing it neither advances nor debits heat; a separate thermal transport
//! owner must account for radiative cooling and internal radiative exchange.

use super::*;
use crate::conductor::{
    CoatingLayer, CoherentAbsorbingCoating, CoherentDielectricCoating, MAX_COHERENT_STACK_LAYERS,
};
use crate::dielectric::DielectricIor;
use crate::thermal_surface::{CoatingEmissionError, CoatingThermalEmission};
use std::collections::BTreeSet;

/// Preserve original native-state and optical/emission refusals.
#[derive(Clone, Debug, PartialEq)]
pub enum NativeCoatingEmissionError {
    Binding(SurfaceHeatingError),
    Emission(CoatingEmissionError),
}
impl From<SurfaceHeatingError> for NativeCoatingEmissionError {
    fn from(value: SurfaceHeatingError) -> Self {
        Self::Binding(value)
    }
}
impl From<DepositionError> for NativeCoatingEmissionError {
    fn from(value: DepositionError) -> Self {
        Self::Binding(value.into())
    }
}
impl From<CoatingEmissionError> for NativeCoatingEmissionError {
    fn from(value: CoatingEmissionError) -> Self {
        Self::Emission(value)
    }
}
impl core::fmt::Display for NativeCoatingEmissionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native coating emission refused: {self:?}")
    }
}
impl std::error::Error for NativeCoatingEmissionError {}

/// Immutable physical observation at one accepted native thermal time.
///
/// Each film uses its own cell temperature for Planck emission and its own
/// existing fixed, temperature or exact-caloric-enthalpy optical law. Distinct
/// enthalpies on a latent plateau therefore retain their optical differences.
/// The low-level owner supplies directional, polarized emission into vacuum/
/// exterior air within the original optical tables' common wavelength domain.
/// No diffuse approximation, RGB lift or bolometric continuation is introduced.
#[derive(Clone, Debug)]
pub struct NativeCoatingEmission {
    emitter: CoatingThermalEmission,
    time_s: f64,
    identity: ContentHash,
}

#[derive(Clone, Copy)]
enum ResolvedFilm {
    Absorbing(CoherentAbsorbingCoating),
    Lossless(CoherentDielectricCoating),
}
impl ResolvedFilm {
    fn layer(&self) -> CoatingLayer<'_> {
        match self {
            Self::Absorbing(optics) => CoatingLayer::Absorbing(optics),
            Self::Lossless(optics) => CoatingLayer::Dielectric(optics),
        }
    }
}

struct Resolver {
    temperatures: BTreeMap<CellId, f64>,
    films: BTreeMap<(CellId, MaterialId, [u8; 32]), CoherentAbsorbingCoating>,
    laws: BTreeSet<(u8, [u8; 32])>,
    rows: usize,
    maximum_rows: usize,
    hash: DomainHasher,
}
impl Resolver {
    fn new(time_s: f64, maximum_rows: usize, source: ContentHash) -> Self {
        let mut hash = DomainHasher::new("org.frankensim.render.native-coating-emission.v1");
        hash.update(source.as_bytes());
        hash.update(&time_s.to_bits().to_le_bytes());
        Self {
            temperatures: BTreeMap::new(),
            films: BTreeMap::new(),
            laws: BTreeSet::new(),
            rows: 0,
            maximum_rows,
            hash,
        }
    }

    fn receiver(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        receiver: SurfaceHeatingReceiver,
    ) -> Result<f64, SurfaceHeatingError> {
        let (lo, hi) = receiver.temperature_bounds_k;
        if !lo.is_finite() || !hi.is_finite() || lo <= 0.0 || hi < lo {
            return Err(invalid(
                "native coating emission needs positive finite temperature support",
            ));
        }
        let cell = thermal
            .cells()
            .get(&receiver.thermal_cell)
            .ok_or_else(|| invalid("native coating emission receiver must be active"))?;
        if cell.material_id != receiver.material_id {
            return Err(invalid(
                "native coating emission must retain its declared material",
            ));
        }
        let temperature = if let Some(&temperature) = self.temperatures.get(&receiver.thermal_cell)
        {
            temperature
        } else {
            let temperature = thermal.temperature_k(receiver.thermal_cell)?;
            self.temperatures.insert(receiver.thermal_cell, temperature);
            temperature
        };
        if !temperature.is_finite() || temperature < lo || temperature > hi {
            return Err(SurfaceHeatingError::TemperatureDomain {
                thermal_cell: receiver.thermal_cell,
                temperature_k: temperature,
                minimum_k: lo,
                maximum_k: hi,
            });
        }
        for id in [receiver.thermal_cell.0, receiver.material_id.0] {
            self.hash.update(&id.to_le_bytes());
        }
        for value in [lo, hi, temperature, cell.specific_h_j_kg] {
            self.hash.update(&value.to_bits().to_le_bytes());
        }
        Ok(temperature)
    }

    fn law(
        &mut self,
        kind: u8,
        identity: ContentHash,
        rows: usize,
    ) -> Result<(), SurfaceHeatingError> {
        self.hash.update(&[kind]);
        self.hash.update(identity.as_bytes());
        if self.laws.insert((kind, *identity.as_bytes())) {
            self.rows = self
                .rows
                .checked_add(rows)
                .ok_or(SurfaceHeatingError::Budget)?;
            if self.rows > self.maximum_rows {
                return Err(SurfaceHeatingError::Budget);
            }
        }
        Ok(())
    }

    fn films(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        films: &[StatefulCoatingLayer],
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Vec<ResolvedFilm>, Vec<f64>), SurfaceHeatingError> {
        if films.len() > MAX_COHERENT_STACK_LAYERS {
            return Err(SurfaceHeatingError::Budget);
        }
        self.hash.update(&(films.len() as u64).to_le_bytes());
        let mut resolved = Vec::with_capacity(films.len());
        let mut temperatures = Vec::with_capacity(films.len());
        for film in films {
            checkpoint()?;
            temperatures.push(self.receiver(thermal, film.receiver)?);
            resolved.push(match &film.optics {
                CoatingOpticsModel::Fixed(optics) => {
                    self.hash.update(&[2]);
                    ResolvedFilm::Absorbing(*optics)
                }
                CoatingOpticsModel::Lossless(optics) => {
                    self.hash.update(&[3]);
                    ResolvedFilm::Lossless(*optics)
                }
                CoatingOpticsModel::Law(law) => {
                    self.law(0, law.identity(), law.row_count())?;
                    let key = (
                        film.receiver.thermal_cell,
                        film.receiver.material_id,
                        *law.identity().as_bytes(),
                    );
                    let optics = if let Some(&optics) = self.films.get(&key) {
                        optics
                    } else {
                        let optics = law.resolve(
                            thermal,
                            film.receiver.thermal_cell,
                            film.receiver.material_id,
                        )?;
                        self.films.insert(key, optics);
                        optics
                    };
                    ResolvedFilm::Absorbing(optics)
                }
            });
        }
        Ok((resolved, temperatures))
    }
}

impl NativeCoatingEmission {
    /// Prepare independently emitting films and their opaque substrate from
    /// original native states. The explicit films replace any substrate coating.
    /// Unique optical-law rows or original material property queries share one
    /// declared preparation budget; repeated
    /// receiver/law pairs resolve once. Fixed optical models need no row budget.
    #[allow(clippy::too_many_arguments)]
    pub fn try_opaque(
        thermal: &DepositionSimulation<'_>,
        films: &[StatefulCoatingLayer],
        substrate_receiver: SurfaceHeatingReceiver,
        substrate: &SubstrateOpticsModel,
        source_identity: ContentHash,
        max_total_optical_knots: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, NativeCoatingEmissionError> {
        checkpoint()?;
        let mut resolver =
            Resolver::new(thermal.time_s(), max_total_optical_knots, source_identity);
        resolver.hash.update(&[0]);
        let (films, temperatures) = resolver.films(thermal, films, &mut checkpoint)?;
        checkpoint()?;
        let substrate_temperature = resolver.receiver(thermal, substrate_receiver)?;
        let substrate = match substrate {
            SubstrateOpticsModel::Fixed(optics) => {
                resolver.hash.update(&[2]);
                *optics
            }
            SubstrateOpticsModel::Law(law) => {
                resolver.law(1, law.identity(), law.budget_size())?;
                law.resolve(
                    thermal,
                    substrate_receiver.thermal_cell,
                    substrate_receiver.material_id,
                )?
            }
        };
        let layers: Vec<_> = films.iter().map(ResolvedFilm::layer).collect();
        let emitter = CoatingThermalEmission::try_opaque(
            &layers,
            &temperatures,
            substrate,
            substrate_temperature,
            source_identity,
        )?;
        Self::finish(thermal.time_s(), emitter, resolver, &mut checkpoint)
    }

    /// Prepare emitting finite films with a lossless opposite medium. The
    /// observed port is vacuum/exterior air; `opposite_ior` names the other
    /// physical port for this view. Reverse traversal retains authored film
    /// temperatures and receivers. A bare interface needs no thermal cells.
    #[allow(clippy::too_many_arguments)]
    pub fn try_transmitting(
        thermal: &DepositionSimulation<'_>,
        films: &[StatefulCoatingLayer],
        opposite_ior: DielectricIor,
        reverse: bool,
        source_identity: ContentHash,
        max_total_optical_knots: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, NativeCoatingEmissionError> {
        checkpoint()?;
        let mut resolver =
            Resolver::new(thermal.time_s(), max_total_optical_knots, source_identity);
        resolver.hash.update(&[1]);
        let (films, temperatures) = resolver.films(thermal, films, &mut checkpoint)?;
        let layers: Vec<_> = films.iter().map(ResolvedFilm::layer).collect();
        let emitter = CoatingThermalEmission::try_transmitting(
            &layers,
            &temperatures,
            opposite_ior,
            reverse,
            source_identity,
        )?;
        Self::finish(thermal.time_s(), emitter, resolver, &mut checkpoint)
    }

    fn finish(
        time_s: f64,
        emitter: CoatingThermalEmission,
        mut resolver: Resolver,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, NativeCoatingEmissionError> {
        checkpoint()?;
        resolver.hash.update(emitter.identity().as_bytes());
        Ok(Self {
            emitter,
            time_s,
            identity: resolver.hash.finalize(),
        })
    }

    /// Reusable native directional/spectral emission owner, retaining its
    /// positive per-film absorption and original polarization convention.
    pub fn emitter(&self) -> &CoatingThermalEmission {
        &self.emitter
    }
    /// Accepted native physical time at which all optical and caloric states agree.
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Exact state, receiver/law bindings and immutable optical-source identity.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}
