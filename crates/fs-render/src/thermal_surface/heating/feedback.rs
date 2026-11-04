//! Source-authored surface optical constants on the original thermal state.
//!
//! The native caloric law alone owns temperature and latent inventory. These
//! tables interpolate absolute substrate n/k, then reconstruct the original
//! conductor owner; no Fresnel intensity, phase fraction or film law is mixed.

use super::SurfaceHeatingError;
use crate::conductor::{ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{CellId, DepositionSimulation, MaterialId, ThermalStorageBasis};
use std::sync::Arc;

#[cfg(feature = "material-state")]
mod material;

#[derive(Clone, Copy, Debug)]
enum Coordinate {
    Temperature,
    SpecificEnthalpy(ContentHash),
}

#[derive(Clone, Debug)]
enum Model {
    Rows {
        source: ConductorSource,
        rows: Arc<[(f64, ConductorOptics)]>,
    },
    #[cfg(feature = "material-state")]
    Material(Arc<material::MaterialTemperatureConductor>),
}

/// Immutable constitutive n/k law versus temperature or specific enthalpy.
///
/// An authored table retains each original nine-knot spectrum and its source.
/// All rows share one spectral grid and exactly the same attached coating.
/// Linear interpolation in the declared state coordinate is an authored law,
/// not an inferred solid/liquid mixing rule or a claim of measured behavior.
/// The reconstructed substrate carries the separately declared law source;
/// its fixed coating keeps its own original source and optical parameters.
/// Wavelength evaluation, coherent films and power balance remain owned by
/// [`ConductorOptics`]. Neither state nor wavelength support is extrapolated.
/// With `material-state`, an original material-card model instead refreshes its
/// pinned claims at the native temperature, without fitting a new state table.
#[derive(Clone, Debug)]
pub struct SurfaceConductorLaw {
    identity: ContentHash,
    coordinate: Coordinate,
    temperature_bounds_k: (f64, f64),
    model: Model,
}

impl SurfaceConductorLaw {
    /// Declare n/k linear in positive kelvin between at least two ordered rows.
    /// A temperature law cannot distinguish latent states at the same T; use
    /// [`Self::try_enthalpy`] when the authored optics resolve those states.
    pub fn try_temperature(
        source: ConductorSource,
        rows: &[(f64, ConductorOptics)],
        max_rows: usize,
    ) -> Result<Self, SurfaceHeatingError> {
        let bounds = rows
            .first()
            .zip(rows.last())
            .map(|(first, last)| (first.0, last.0))
            .ok_or_else(|| invalid("surface optical law needs at least two rows"))?;
        Self::admit(source, Coordinate::Temperature, bounds, rows, max_rows)
    }

    /// Declare n/k linear in the native caloric law's specific enthalpy [J/kg].
    /// Negative coordinates are valid. The nonzero caloric identity binds the
    /// exact law and energy datum; equal-temperature latent states stay distinct.
    /// The independent positive-kelvin domain remains binding, including on a
    /// plateau where its lower and upper temperature may be equal.
    pub fn try_enthalpy(
        source: ConductorSource,
        caloric_identity: ContentHash,
        temperature_bounds_k: (f64, f64),
        rows: &[(f64, ConductorOptics)],
        max_rows: usize,
    ) -> Result<Self, SurfaceHeatingError> {
        if caloric_identity == ContentHash([0; 32]) {
            return Err(invalid(
                "surface enthalpy optics require a caloric identity",
            ));
        }
        Self::admit(
            source,
            Coordinate::SpecificEnthalpy(caloric_identity),
            temperature_bounds_k,
            rows,
            max_rows,
        )
    }

    fn admit(
        source: ConductorSource,
        coordinate: Coordinate,
        temperature_bounds_k: (f64, f64),
        rows: &[(f64, ConductorOptics)],
        max_rows: usize,
    ) -> Result<Self, SurfaceHeatingError> {
        if rows.len() > max_rows {
            return Err(SurfaceHeatingError::Budget);
        }
        let (minimum_k, maximum_k) = temperature_bounds_k;
        if rows.len() < 2
            || !minimum_k.is_finite()
            || !maximum_k.is_finite()
            || minimum_k <= 0.0
            || maximum_k < minimum_k
            || rows.iter().any(|(x, _)| !x.is_finite())
            || rows.windows(2).any(|pair| {
                let span = pair[1].0 - pair[0].0;
                !span.is_finite() || span <= 0.0
            })
        {
            return Err(invalid(
                "ordered finite optical rows and positive kelvin support required",
            ));
        }
        let reference = rows[0].1;
        let wavelengths = reference
            .samples()
            .map(|sample| sample.wavelength_nm().to_bits());
        for (_, optics) in rows {
            if optics
                .samples()
                .map(|sample| sample.wavelength_nm().to_bits())
                != wavelengths
                || optics.max_knot_span_nm().to_bits() != reference.max_knot_span_nm().to_bits()
                || optics.coating() != reference.coating()
                || optics.absorbing_coating() != reference.absorbing_coating()
                || optics.periodic_coating() != reference.periodic_coating()
            {
                return Err(invalid(
                    "surface optical rows must share their spectral grid, gap limit and fixed coating",
                ));
            }
            if source.status() == ConductorDataStatus::CallerAssertedMeasured
                && optics.source().status() != ConductorDataStatus::CallerAssertedMeasured
            {
                return Err(invalid(
                    "surface interpolation cannot upgrade endpoint data to measured optics",
                ));
            }
        }
        let mut hash = DomainHasher::new("org.frankensim.render.surface-conductor-law.v1");
        hash.update(source.identity().as_bytes());
        hash.update(&[match source.status() {
            ConductorDataStatus::Representative => 0,
            ConductorDataStatus::CallerAssertedMeasured => 1,
            ConductorDataStatus::MaterialStateClaims => 2,
        }]);
        match coordinate {
            Coordinate::Temperature => hash.update(&[0]),
            Coordinate::SpecificEnthalpy(caloric) => {
                hash.update(&[1]);
                hash.update(caloric.as_bytes());
            }
        }
        for value in [minimum_k, maximum_k] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.update(&(rows.len() as u64).to_le_bytes());
        for (x, optics) in rows {
            hash.update(&x.to_bits().to_le_bytes());
            hash.update(optics.content_identity().as_bytes());
        }
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(rows.len())
            .map_err(|_| SurfaceHeatingError::Budget)?;
        retained.extend_from_slice(rows);
        Ok(Self {
            identity: hash.finalize(),
            coordinate,
            temperature_bounds_k,
            model: Model::Rows {
                source,
                rows: retained.into(),
            },
        })
    }

    /// Exact constitutive rows, source claims, coordinate and support identity.
    pub fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Authored state rows. Original material models have no fitted rows.
    pub fn row_count(&self) -> usize {
        match &self.model {
            Model::Rows { rows, .. } => rows.len(),
            #[cfg(feature = "material-state")]
            Model::Material(_) => 0,
        }
    }

    /// Storage/query budget units: authored rows or all retained source properties.
    pub(crate) fn budget_size(&self) -> usize {
        match &self.model {
            Model::Rows { rows, .. } => rows.len(),
            #[cfg(feature = "material-state")]
            Model::Material(source) => source.property_count(),
        }
    }

    pub fn temperature_bounds_k(&self) -> (f64, f64) {
        self.temperature_bounds_k
    }

    pub fn caloric_identity(&self) -> Option<ContentHash> {
        match self.coordinate {
            Coordinate::Temperature => None,
            Coordinate::SpecificEnthalpy(identity) => Some(identity),
        }
    }

    pub fn specific_enthalpy_bounds_j_kg(&self) -> Option<(f64, f64)> {
        self.caloric_identity().map(|_| self.coordinate_bounds())
    }

    /// Evaluate at an active native cell, retaining its actual material and,
    /// for an h law, the exact original storage basis, caloric identity and datum.
    /// Native temperature and h are checked even if the surface is unilluminated.
    pub fn resolve(
        &self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        material: MaterialId,
    ) -> Result<ConductorOptics, SurfaceHeatingError> {
        self.optics_at_coordinate(self.coordinate_at(thermal, cell, material)?)
    }

    /// Validate the complete native binding without reconstructing n/k. Used
    /// by accepted-endpoint checks and shared-cell response caches.
    pub(crate) fn coordinate_at(
        &self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        material: MaterialId,
    ) -> Result<f64, SurfaceHeatingError> {
        let native = thermal
            .cells()
            .get(&cell)
            .ok_or_else(|| invalid("surface optical receiver must be an active native cell"))?;
        if native.material_id != material {
            return Err(invalid("surface optical receiver material changed"));
        }
        let temperature_k = thermal.temperature_k(cell)?;
        if !self.supports_temperature(temperature_k) {
            return Err(SurfaceHeatingError::TemperatureDomain {
                thermal_cell: cell,
                temperature_k,
                minimum_k: self.temperature_bounds_k.0,
                maximum_k: self.temperature_bounds_k.1,
            });
        }
        let coordinate = match self.coordinate {
            Coordinate::Temperature => temperature_k,
            Coordinate::SpecificEnthalpy(identity) => {
                if thermal.material_storage_basis(material)? != ThermalStorageBasis::Enthalpy
                    || thermal.material_caloric_identity(material)? != Some(*identity.as_bytes())
                {
                    return Err(invalid(
                        "surface enthalpy optics require the original caloric identity and enthalpy storage",
                    ));
                }
                native.specific_h_j_kg
            }
        };
        self.check_coordinate(coordinate)?;
        Ok(coordinate)
    }

    /// Endpoint validation includes the original source queries, even in a
    /// dark interval where no new Fresnel packet is needed. Row models need
    /// only their already-admitted coordinate support checked. A material model
    /// returns its freshly admitted optics so each bound incidence can validate
    /// its wavelengths without querying the source or solving Fresnel again.
    pub(crate) fn validate_native(
        &self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        material: MaterialId,
    ) -> Result<Option<ConductorOptics>, SurfaceHeatingError> {
        let coordinate = self.coordinate_at(thermal, cell, material)?;
        #[cfg(feature = "material-state")]
        if let Model::Material(source) = &self.model {
            return source.optics_at_temperature(coordinate).map(Some);
        }
        #[cfg(not(feature = "material-state"))]
        let _ = coordinate;
        Ok(None)
    }

    /// Evaluate a temperature-coordinate law directly. An h-coordinate law
    /// refuses because temperature does not identify a state within latent heat.
    /// Native heating callers should use [`Self::resolve`] for material binding.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ConductorOptics, SurfaceHeatingError> {
        if !matches!(self.coordinate, Coordinate::Temperature)
            || !self.supports_temperature(temperature_k)
        {
            return Err(invalid(
                "surface law requires its declared temperature coordinate and support",
            ));
        }
        self.optics_at_coordinate(temperature_k)
    }

    /// Evaluate an h-coordinate law directly in its declared caloric datum.
    /// The caller supplies temperature from that same caloric law; use
    /// [`Self::resolve`] to check the native material, basis and caloric identity.
    /// Temperature-coordinate laws refuse this physically different coordinate.
    pub fn at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        temperature_k: f64,
    ) -> Result<ConductorOptics, SurfaceHeatingError> {
        if !matches!(self.coordinate, Coordinate::SpecificEnthalpy(_))
            || !self.supports_temperature(temperature_k)
        {
            return Err(invalid(
                "surface law requires its declared enthalpy coordinate and temperature support",
            ));
        }
        self.optics_at_coordinate(specific_h_j_kg)
    }

    fn supports_temperature(&self, temperature_k: f64) -> bool {
        temperature_k.is_finite()
            && (self.temperature_bounds_k.0..=self.temperature_bounds_k.1).contains(&temperature_k)
    }

    fn check_coordinate(&self, coordinate: f64) -> Result<(), SurfaceHeatingError> {
        let (minimum, maximum) = self.coordinate_bounds();
        if !coordinate.is_finite() || !(minimum..=maximum).contains(&coordinate) {
            return Err(invalid(
                "surface optical state is outside its constitutive coordinate support",
            ));
        }
        Ok(())
    }

    fn coordinate_bounds(&self) -> (f64, f64) {
        match &self.model {
            Model::Rows { rows, .. } => (rows[0].0, rows[rows.len() - 1].0),
            #[cfg(feature = "material-state")]
            Model::Material(_) => self.temperature_bounds_k,
        }
    }

    pub(crate) fn optics_at_coordinate(
        &self,
        coordinate: f64,
    ) -> Result<ConductorOptics, SurfaceHeatingError> {
        self.check_coordinate(coordinate)?;
        let (source, rows) = match &self.model {
            Model::Rows { source, rows } => (*source, rows),
            #[cfg(feature = "material-state")]
            Model::Material(source) => return source.optics_at_temperature(coordinate),
        };
        let upper = rows
            .partition_point(|(x, _)| *x < coordinate)
            .clamp(1, rows.len() - 1);
        let (a, left) = rows[upper - 1];
        let (b, right) = rows[upper];
        let fraction = (coordinate - a) / (b - a);
        let mut samples = left.samples();
        for (sample, rhs) in samples.iter_mut().zip(right.samples()) {
            *sample = ConductorIorSample::try_new_in_band(
                sample.wavelength_nm(),
                interpolate_value(sample.eta(), rhs.eta(), fraction),
                interpolate_value(sample.k(), rhs.k(), fraction),
            )?;
        }
        let mut optics =
            ConductorOptics::try_new_in_band(samples, left.max_knot_span_nm(), source)?;
        if let Some(coating) = left.coating() {
            optics = optics.with_coating(coating);
        } else if let Some(coating) = left.absorbing_coating() {
            optics = optics.with_absorbing_coating(coating);
        } else if let Some(coating) = left.periodic_coating() {
            optics = optics.with_periodic_coating(coating);
        }
        Ok(optics)
    }
}

fn interpolate_value(a: f64, b: f64, fraction: f64) -> f64 {
    if a == b || fraction == 0.0 {
        a
    } else if fraction == 1.0 {
        b
    } else if fraction <= 0.5 {
        // n/k are positive and at most 100, so subtraction cannot overflow.
        // One rounding preserves subnormal contributions lost by separately
        // rounded convex products; the nearer endpoint limits cancellation.
        fraction.mul_add(b - a, a)
    } else {
        (1.0 - fraction).mul_add(a - b, b)
    }
}

fn invalid(message: &'static str) -> SurfaceHeatingError {
    SurfaceHeatingError::Invalid(message)
}

#[cfg(test)]
mod tests;
