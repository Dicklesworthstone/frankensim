//! State-dependent opaque diffuse absorptivity on the original thermal state.
//!
//! Kirchhoff reciprocity is retained explicitly: the same authored spectral
//! emissivity is the surface absorptance used by incident heating. The native
//! caloric law alone owns temperature and latent inventory. Tables interpolate
//! emissivity in the declared state coordinate before wavelength interpolation;
//! no RGB lift, phase mixing rule or infrared continuation is invented.

use super::{SurfaceHeatingError, SurfacePowerResponse};
use crate::thermal_surface::{SpectralDiffuseThermalSurface, SpectralEmissivity};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{CellId, DepositionSimulation, MaterialId, ThermalStorageBasis};
use std::sync::Arc;

#[derive(Clone, Copy, Debug)]
enum Coordinate {
    Temperature,
    SpecificEnthalpy(ContentHash),
}

/// Immutable spectral emissivity versus native temperature or specific enthalpy.
///
/// Every row must use the same visible wavelength knots. Linear interpolation in
/// state is an authored constitutive law. Wavelength interpolation remains owned
/// by [`SpectralEmissivity`]. Enthalpy laws bind the exact caloric identity and
/// datum, so equal-temperature states inside latent heat may have different
/// absorptivity. Neither state nor wavelength support is extrapolated.
#[derive(Clone, Debug)]
pub struct SurfaceDiffuseLaw {
    identity: ContentHash,
    coordinate: Coordinate,
    temperature_bounds_k: (f64, f64),
    rows: Arc<[(f64, SpectralEmissivity)]>,
}

impl SurfaceDiffuseLaw {
    pub fn try_temperature(
        source_identity: ContentHash,
        rows: &[(f64, SpectralEmissivity)],
        max_rows: usize,
    ) -> Result<Self, SurfaceHeatingError> {
        let bounds = rows
            .first()
            .zip(rows.last())
            .map(|(first, last)| (first.0, last.0))
            .ok_or_else(|| invalid("diffuse optical law needs at least two rows"))?;
        Self::admit(
            source_identity,
            Coordinate::Temperature,
            bounds,
            rows,
            max_rows,
        )
    }

    pub fn try_enthalpy(
        source_identity: ContentHash,
        caloric_identity: ContentHash,
        temperature_bounds_k: (f64, f64),
        rows: &[(f64, SpectralEmissivity)],
        max_rows: usize,
    ) -> Result<Self, SurfaceHeatingError> {
        if caloric_identity == ContentHash([0; 32]) {
            return Err(invalid(
                "diffuse enthalpy optics require a caloric identity",
            ));
        }
        Self::admit(
            source_identity,
            Coordinate::SpecificEnthalpy(caloric_identity),
            temperature_bounds_k,
            rows,
            max_rows,
        )
    }

    fn admit(
        source_identity: ContentHash,
        coordinate: Coordinate,
        temperature_bounds_k: (f64, f64),
        rows: &[(f64, SpectralEmissivity)],
        max_rows: usize,
    ) -> Result<Self, SurfaceHeatingError> {
        if source_identity == ContentHash([0; 32]) {
            return Err(invalid("diffuse optical law requires a source identity"));
        }
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
                "ordered diffuse optical rows and positive kelvin support required",
            ));
        }
        let reference = rows[0].1.knots();
        if rows.iter().skip(1).any(|(_, emissivity)| {
            emissivity.knots().len() != reference.len()
                || emissivity
                    .knots()
                    .iter()
                    .zip(reference)
                    .any(|(a, b)| a.0.to_bits() != b.0.to_bits())
        }) {
            return Err(invalid(
                "diffuse optical rows must share one spectral wavelength grid",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.render.surface-diffuse-law.v1");
        hash.update(source_identity.as_bytes());
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
        for (x, emissivity) in rows {
            hash.update(&x.to_bits().to_le_bytes());
            hash.update(emissivity.identity().as_bytes());
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
            rows: retained.into(),
        })
    }

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

    pub fn resolve(
        &self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        material: MaterialId,
    ) -> Result<SpectralEmissivity, SurfaceHeatingError> {
        self.emissivity_at_coordinate(self.coordinate_at(thermal, cell, material)?)
    }

    pub fn resolve_response<const N: usize>(
        &self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        material: MaterialId,
        wavelengths_nm: [f64; N],
    ) -> Result<SurfacePowerResponse<N>, SurfaceHeatingError> {
        let emissivity = self.resolve(thermal, cell, material)?;
        SurfacePowerResponse::from_spectral_diffuse(
            wavelengths_nm,
            &emissivity,
            self.identity,
        )
    }

    /// Resolve the reciprocal equilibrium-emission view from the SAME native
    /// state law used by incident absorption. The native caloric owner supplies
    /// temperature; no separate emissivity fit or phase-coordinate mapping is
    /// introduced. This remains a visible-band source and makes no bolometric
    /// cooling claim.
    pub fn resolve_surface(
        &self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        material: MaterialId,
    ) -> Result<SpectralDiffuseThermalSurface, SurfaceHeatingError> {
        let emissivity = self.resolve(thermal, cell, material)?;
        Ok(SpectralDiffuseThermalSurface::try_new(
            thermal.temperature_k(cell)?,
            emissivity,
            self.identity,
        )?)
    }

    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<SpectralEmissivity, SurfaceHeatingError> {
        if !matches!(self.coordinate, Coordinate::Temperature)
            || !self.supports_temperature(temperature_k)
        {
            return Err(invalid(
                "diffuse law requires its declared temperature coordinate and support",
            ));
        }
        self.emissivity_at_coordinate(temperature_k)
    }

    pub fn at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        temperature_k: f64,
    ) -> Result<SpectralEmissivity, SurfaceHeatingError> {
        if !matches!(self.coordinate, Coordinate::SpecificEnthalpy(_))
            || !self.supports_temperature(temperature_k)
        {
            return Err(invalid(
                "diffuse law requires its declared enthalpy coordinate and temperature support",
            ));
        }
        self.emissivity_at_coordinate(specific_h_j_kg)
    }

    pub(crate) fn coordinate_at(
        &self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        material: MaterialId,
    ) -> Result<f64, SurfaceHeatingError> {
        let native = thermal
            .cells()
            .get(&cell)
            .ok_or_else(|| invalid("diffuse optical receiver must be an active native cell"))?;
        if native.material_id != material {
            return Err(invalid("diffuse optical receiver material changed"));
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
                        "diffuse enthalpy optics require the original caloric identity and enthalpy storage",
                    ));
                }
                native.specific_h_j_kg
            }
        };
        self.check_coordinate(coordinate)?;
        Ok(coordinate)
    }

    fn supports_temperature(&self, temperature_k: f64) -> bool {
        temperature_k.is_finite()
            && (self.temperature_bounds_k.0..=self.temperature_bounds_k.1)
                .contains(&temperature_k)
    }

    fn check_coordinate(&self, coordinate: f64) -> Result<(), SurfaceHeatingError> {
        if !coordinate.is_finite()
            || !(self.rows[0].0..=self.rows[self.rows.len() - 1].0).contains(&coordinate)
        {
            return Err(invalid(
                "diffuse optical state is outside its constitutive coordinate support",
            ));
        }
        Ok(())
    }

    pub(crate) fn emissivity_at_coordinate(
        &self,
        coordinate: f64,
    ) -> Result<SpectralEmissivity, SurfaceHeatingError> {
        self.check_coordinate(coordinate)?;
        let upper = self
            .rows
            .partition_point(|(x, _)| *x < coordinate)
            .clamp(1, self.rows.len() - 1);
        let (a, left) = self.rows[upper - 1];
        let (b, right) = self.rows[upper];
        let fraction = (coordinate - a) / (b - a);
        let mut knots = Vec::new();
        knots
            .try_reserve_exact(left.knots().len())
            .map_err(|_| SurfaceHeatingError::Budget)?;
        for (&(wavelength_nm, emissivity), &(_, rhs)) in
            left.knots().iter().zip(right.knots())
        {
            knots.push((
                wavelength_nm,
                interpolate_value(emissivity, rhs, fraction),
            ));
        }
        SpectralEmissivity::table(&knots).map_err(Into::into)
    }
}

fn interpolate_value(a: f64, b: f64, fraction: f64) -> f64 {
    if a == b || fraction == 0.0 {
        a
    } else if fraction == 1.0 {
        b
    } else if fraction <= 0.5 {
        fraction.mul_add(b - a, a)
    } else {
        (1.0 - fraction).mul_add(a - b, b)
    }
}

fn invalid(message: &'static str) -> SurfaceHeatingError {
    SurfaceHeatingError::Invalid(message)
}
