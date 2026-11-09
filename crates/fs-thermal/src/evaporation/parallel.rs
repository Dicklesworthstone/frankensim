//! Heterogeneous exposed films sharing one original liquid inventory.
//!
//! Each patch sees the same liquid temperature and far-field gas pressures and
//! film temperature, but retains its own exposed area and molecular transfer
//! coefficient. The area-mean flux is `sum(A_i j_i) / sum(A_i)`. Every `j_i`
//! closes its own original nonlinear Stefan film and selected interface kinetics.
//! Averaging the coefficients first is generally wrong for finite kinetics.
//!
//! This wrapper owns no mass, caloric chart or gas. The enclosing owner uses
//! [`ParallelFilmEvaporation::area_m2`] for its ONE liquid inventory, then stores
//! the resulting mass and heat once. Patches are immutable and borrowed; work
//! is bounded by their slice, with no allocation during construction or flux
//! evaluation. Uniform coefficients take the original single-film fast path.

use super::{
    EvaporationError, EvaporationLaw, EvaporationProperties, GasFilm, InterfaceKinetics,
    SurfaceFlux, checked, nonnegative, positive,
};
use crate::enthalpy::{
    DepositedCell, DepositionError, ThermalLaw, ThermalStorageBasis,
    advection::species::{CaloricComposition, CaloricEntropyAnchor},
};

/// One actually exposed patch of the SAME pure-liquid inventory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvaporationPatch {
    pub area_m2: f64,
    /// Molecular coefficient before the Stefan correction; zero seals this patch.
    pub mass_transfer_coefficient_m_s: f64,
}

/// Borrow the original source and retain separate resistances on its patches.
///
/// All patches share one liquid temperature, phase law, kinetic selection and
/// far-field gas. Each patch replaces the incoming [`GasFilm`] coefficient;
/// incoming total/volatile pressures and representative temperature are retained.
/// Geometry/connectivity and distinct face identities belong to the caller.
/// The caller must use [`Self::area_m2`] as the exposed area of its one inventory.
pub struct ParallelFilmEvaporation<'a, L: EvaporationLaw + ?Sized> {
    source: &'a L,
    patches: &'a [EvaporationPatch],
    area_m2: f64,
    uniform_coefficient_m_s: Option<f64>,
}
impl<L: EvaporationLaw + ?Sized> Copy for ParallelFilmEvaporation<'_, L> {}
impl<L: EvaporationLaw + ?Sized> Clone for ParallelFilmEvaporation<'_, L> {
    fn clone(&self) -> Self {
        *self
    }
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), EvaporationError> {
        let (sum, error) = fs_math::eft::two_sum(self.value, value);
        self.value = checked(sum)?;
        self.correction = checked(self.correction + error)?;
        Ok(())
    }
    fn total(&self) -> Result<f64, EvaporationError> {
        checked(self.value + self.correction)
    }
}

impl<'a, L: EvaporationLaw + ?Sized> ParallelFilmEvaporation<'a, L> {
    /// Admit finite nonnegative patch areas/coefficients and positive total area.
    /// Zero-area patches retain validation but contribute no flux or source call.
    pub fn try_new(
        source: &'a L,
        patches: &'a [EvaporationPatch],
    ) -> Result<Self, EvaporationError> {
        let mut total = Sum::default();
        let mut first_coefficient = None;
        let mut uniform = true;
        for patch in patches {
            if !nonnegative(patch.area_m2) || !nonnegative(patch.mass_transfer_coefficient_m_s) {
                return Err(EvaporationError::Invalid(
                    "parallel films require finite nonnegative areas and coefficients",
                ));
            }
            total.add(patch.area_m2)?;
            if patch.area_m2 > 0.0 {
                let first = *first_coefficient.get_or_insert(patch.mass_transfer_coefficient_m_s);
                uniform &= first == patch.mass_transfer_coefficient_m_s;
            }
        }
        let area_m2 = total.total()?;
        if !positive(area_m2) {
            return Err(EvaporationError::Invalid(
                "parallel films require positive total exposed area",
            ));
        }
        Ok(Self {
            source,
            patches,
            area_m2,
            uniform_coefficient_m_s: if uniform { first_coefficient } else { None },
        })
    }

    pub const fn source(&self) -> &'a L {
        self.source
    }
    pub const fn patches(&self) -> &'a [EvaporationPatch] {
        self.patches
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
}

// Normalize before summation to avoid overflowing A*j when the mean is finite.
// Alternative orders recover a tiny area fraction when its final contribution
// is representable even though A/A_total alone underflows.
fn weighted(value: f64, area: f64, total: f64) -> Result<f64, EvaporationError> {
    let direct = (area / total) * value;
    if value == 0.0 || direct != 0.0 {
        return checked(direct);
    }
    for alternate in [(area * value) / total, (value / total) * area] {
        if alternate.is_finite() && alternate != 0.0 {
            return Ok(alternate);
        }
    }
    Err(EvaporationError::Unrepresentable(
        "nonzero parallel-film contribution underflows",
    ))
}

impl<L: EvaporationLaw + ?Sized> EvaporationLaw for ParallelFilmEvaporation<'_, L> {
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
        film: GasFilm,
        kinetics: InterfaceKinetics,
    ) -> Result<SurfaceFlux, EvaporationError> {
        film.validate()?;
        if let Some(coefficient) = self.uniform_coefficient_m_s {
            return self.source.evaluate_film(
                admitted_properties,
                surface_temperature_k,
                GasFilm {
                    mass_transfer_coefficient_m_s: coefficient,
                    ..film
                },
                kinetics,
            );
        }
        let mut mass = Sum::default();
        let mut derivative = Sum::default();
        for patch in self.patches {
            if patch.area_m2 == 0.0 {
                continue;
            }
            let flux = self.source.evaluate_film(
                admitted_properties,
                surface_temperature_k,
                GasFilm {
                    mass_transfer_coefficient_m_s: patch.mass_transfer_coefficient_m_s,
                    ..film
                },
                kinetics,
            )?;
            mass.add(weighted(
                flux.mass_flux_kg_m2_s,
                patch.area_m2,
                self.area_m2,
            )?)?;
            derivative.add(weighted(
                flux.mass_flux_derivative_kg_m2_s_k,
                patch.area_m2,
                self.area_m2,
            )?)?;
        }
        let mass = mass.total()?;
        // One original phase law applies to every patch. Derive these accounts
        // from the SAME represented mean mass instead of rounding four totals.
        Ok(SurfaceFlux {
            mass_flux_kg_m2_s: mass,
            mass_flux_derivative_kg_m2_s_k: derivative.total()?,
            latent_heat_flux_w_m2: checked(mass * admitted_properties.latent_heat_j_kg)?,
            outward_normal_velocity_m_s: checked(-mass / admitted_properties.liquid_density_kg_m3)?,
        })
    }
}

impl<L: EvaporationLaw + ThermalLaw + ?Sized> ThermalLaw for ParallelFilmEvaporation<'_, L> {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.source.storage_basis()
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.source.caloric_identity()
    }
    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        self.source.validate_cell(cell)
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.source.temperature_k(h)
    }
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        self.source.bulk_density_kg_m3(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.source.maximum_temperature_slope()
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h_j_kg: f64,
        final_h_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.source
            .specific_entropy_change_j_kg_k(initial_h_j_kg, final_h_j_kg)
    }
    fn advected_composition(&self) -> Option<CaloricComposition<'_>> {
        self.source.advected_composition()
    }
    fn advected_entropy_anchor(&self) -> Option<CaloricEntropyAnchor<'_>> {
        self.source.advected_entropy_anchor()
    }
}

#[cfg(test)]
mod tests;
