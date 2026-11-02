//! Source-authored film optical constants on the original thermal state.
//!
//! The native caloric law owns temperature and latent inventory. These tables
//! interpolate absolute n/k before the original coherent optical solve; they
//! do not mix absorbed intensities or infer phase fractions or thermal expansion.

use super::SurfaceHeatingError;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentCoatingStack,
    ConductorDataStatus, ConductorSource,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{CellId, DepositionSimulation, MaterialId, ThermalStorageBasis};
use std::sync::Arc;

#[derive(Clone, Copy, Debug)]
enum Coordinate {
    Temperature,
    SpecificEnthalpy(ContentHash),
}

/// An explicit constitutive law for one uniform passive coherent film.
///
/// All rows retain the same nine-knot spectral grid and physical thickness.
/// Linear n/k interpolation is a caller-authored state law, not an inferred
/// solid/liquid mixing rule or a claim of measured behavior. Exact zero
/// extinction remains transparent. The reconstructed film carries the declared
/// law source; identities retain every endpoint's original source and status.
///
/// [CoherentAbsorbingCoating] continues to own wavelength interpolation,
/// passivity and coherent propagation. Neither state nor wavelength support
/// is extrapolated. Evaluation allocates no memory; clones share admitted rows.
#[derive(Clone, Debug)]
pub struct SurfaceCoatingLaw {
    identity: ContentHash,
    source: ConductorSource,
    coordinate: Coordinate,
    temperature_bounds_k: (f64, f64),
    rows: Arc<[(f64, CoherentAbsorbingCoating)]>,
}

impl SurfaceCoatingLaw {
    /// Declare n/k linear in positive kelvin between at least two ordered rows.
    /// Use [Self::try_enthalpy] when the authored film changes between latent
    /// states at the same temperature.
    pub fn try_temperature(
        source: ConductorSource,
        rows: &[(f64, CoherentAbsorbingCoating)],
        max_rows: usize,
    ) -> Result<Self, SurfaceHeatingError> {
        let bounds = rows
            .first()
            .zip(rows.last())
            .map(|(first, last)| (first.0, last.0))
            .ok_or_else(|| invalid("surface coating law needs at least two rows"))?;
        Self::admit(source, Coordinate::Temperature, bounds, rows, max_rows)
    }

    /// Declare n/k linear in the native caloric law's specific enthalpy [J/kg].
    /// Negative coordinates are valid. The nonzero caloric identity binds the
    /// exact law and its energy datum. The independent positive-kelvin support
    /// remains binding, including a latent plateau with equal temperature bounds.
    pub fn try_enthalpy(
        source: ConductorSource,
        caloric_identity: ContentHash,
        temperature_bounds_k: (f64, f64),
        rows: &[(f64, CoherentAbsorbingCoating)],
        max_rows: usize,
    ) -> Result<Self, SurfaceHeatingError> {
        if caloric_identity == ContentHash([0; 32]) {
            return Err(invalid(
                "surface enthalpy coating requires a caloric identity",
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
        rows: &[(f64, CoherentAbsorbingCoating)],
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
                "ordered finite coating rows and positive kelvin support required",
            ));
        }
        let reference = rows[0].1;
        let wavelengths = reference
            .samples()
            .map(|sample| sample.wavelength_nm().to_bits());
        for (_, coating) in rows {
            if coating
                .samples()
                .map(|sample| sample.wavelength_nm().to_bits())
                != wavelengths
                || coating.thickness_m().to_bits() != reference.thickness_m().to_bits()
            {
                return Err(invalid(
                    "surface coating rows must share their spectral grid and fixed thickness",
                ));
            }
            if source.status() == ConductorDataStatus::CallerAssertedMeasured
                && coating.source().status() != ConductorDataStatus::CallerAssertedMeasured
            {
                return Err(invalid(
                    "surface interpolation cannot upgrade endpoint data to measured optics",
                ));
            }
        }
        let mut hash = DomainHasher::new("org.frankensim.render.surface-coating-law.v1");
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
        for (x, coating) in rows {
            hash.update(&x.to_bits().to_le_bytes());
            // Reuse the original film identity, including its source status.
            let layers = [CoatingLayer::Absorbing(coating)];
            hash.update(
                CoherentCoatingStack::try_new(&layers)?
                    .content_identity()
                    .as_bytes(),
            );
        }
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(rows.len())
            .map_err(|_| SurfaceHeatingError::Budget)?;
        retained.extend_from_slice(rows);
        Ok(Self {
            identity: hash.finalize(),
            source,
            coordinate,
            temperature_bounds_k,
            rows: retained.into(),
        })
    }

    /// Exact constitutive rows, source claims, coordinate and support identity.
    pub fn identity(&self) -> ContentHash {
        self.identity
    }

    pub fn row_count(&self) -> usize {
        self.rows.len()
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
        self.caloric_identity()
            .map(|_| (self.rows[0].0, self.rows[self.rows.len() - 1].0))
    }

    /// Evaluate at an active native cell with its actual material and, for an
    /// h law, the original storage basis, caloric identity and energy datum.
    /// Native temperature and h are checked even when the film is unilluminated.
    pub fn resolve(
        &self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        material: MaterialId,
    ) -> Result<CoherentAbsorbingCoating, SurfaceHeatingError> {
        self.coating_at_coordinate(self.coordinate_at(thermal, cell, material)?)
    }

    /// Validate native binding without reconstructing the coating spectrum.
    pub(crate) fn coordinate_at(
        &self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        material: MaterialId,
    ) -> Result<f64, SurfaceHeatingError> {
        let native = thermal
            .cells()
            .get(&cell)
            .ok_or_else(|| invalid("surface coating receiver must be an active native cell"))?;
        if native.material_id != material {
            return Err(invalid("surface coating receiver material changed"));
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
                        "surface enthalpy coating requires the original caloric identity and enthalpy storage",
                    ));
                }
                native.specific_h_j_kg
            }
        };
        self.check_coordinate(coordinate)?;
        Ok(coordinate)
    }

    /// Evaluate a temperature-coordinate law directly. An h-coordinate law
    /// refuses because temperature alone cannot resolve latent states. Use
    /// [Self::resolve] to check the actual native receiver and material.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<CoherentAbsorbingCoating, SurfaceHeatingError> {
        if !matches!(self.coordinate, Coordinate::Temperature)
            || !self.supports_temperature(temperature_k)
        {
            return Err(invalid(
                "surface coating law requires its declared temperature coordinate and support",
            ));
        }
        self.coating_at_coordinate(temperature_k)
    }

    /// Evaluate an h-coordinate law in its declared caloric datum. The caller
    /// supplies temperature from that same caloric law; [Self::resolve] checks
    /// the actual native material, storage basis and caloric identity.
    pub fn at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        temperature_k: f64,
    ) -> Result<CoherentAbsorbingCoating, SurfaceHeatingError> {
        if !matches!(self.coordinate, Coordinate::SpecificEnthalpy(_))
            || !self.supports_temperature(temperature_k)
        {
            return Err(invalid(
                "surface coating law requires its declared enthalpy coordinate and temperature support",
            ));
        }
        self.coating_at_coordinate(specific_h_j_kg)
    }

    fn supports_temperature(&self, temperature_k: f64) -> bool {
        temperature_k.is_finite()
            && (self.temperature_bounds_k.0..=self.temperature_bounds_k.1).contains(&temperature_k)
    }

    fn check_coordinate(&self, coordinate: f64) -> Result<(), SurfaceHeatingError> {
        if !coordinate.is_finite()
            || !(self.rows[0].0..=self.rows[self.rows.len() - 1].0).contains(&coordinate)
        {
            return Err(invalid(
                "surface coating state is outside its constitutive coordinate support",
            ));
        }
        Ok(())
    }

    pub(crate) fn coating_at_coordinate(
        &self,
        coordinate: f64,
    ) -> Result<CoherentAbsorbingCoating, SurfaceHeatingError> {
        self.check_coordinate(coordinate)?;
        let upper = self
            .rows
            .partition_point(|(x, _)| *x < coordinate)
            .clamp(1, self.rows.len() - 1);
        let (a, left) = self.rows[upper - 1];
        let (b, right) = self.rows[upper];
        let fraction = (coordinate - a) / (b - a);
        let mut samples = left.samples();
        for (sample, rhs) in samples.iter_mut().zip(right.samples()) {
            *sample = CoatingIorSample::try_new(
                sample.wavelength_nm(),
                interpolate_value(sample.eta(), rhs.eta(), fraction),
                interpolate_value(sample.k(), rhs.k(), fraction),
            )?;
        }
        Ok(CoherentAbsorbingCoating::try_new(
            samples,
            left.thickness_m(),
            self.source,
        )?)
    }
}

fn interpolate_value(a: f64, b: f64, fraction: f64) -> f64 {
    if a == b || fraction == 0.0 {
        a
    } else if fraction == 1.0 {
        b
    } else if fraction <= 0.5 {
        // n/k are nonnegative and at most 100. A single rounding retains
        // subnormal contributions; the nearer endpoint limits cancellation.
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
