//! Original dielectric claims on the accepted native slab temperature clock.

use super::{SurfaceHeatingError, invalid};
use crate::material_state::MaterialOpticalBinding;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::MaterialCard;
use fs_material::state_point::{MaterialStatePointError, SpectralOpticalStatePoint};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use std::sync::Arc;

/// A homogeneous slab's original phase-index and intensity-absorption claims.
///
/// Temperature changes are evaluated by the original material claim owner,
/// including every selected auxiliary property, wavelength and fit-domain
/// limit. The original family, selected claims and all non-temperature axes
/// remain fixed. No temperature fit, emissivity fit or extrapolation is added.
/// Tabulated phase index and analytic Sellmeier dispersion retain their native
/// parameterizations and independent absorption spectra.
///
/// This is the weakly absorbing, real-interface incoherent slab model, not a
/// coherent complex-index film. Geometry, opposite-medium index, isothermal
/// storage and thermal/material association are separately declared by its
/// native receiver. No expansion, density update or phase mixture is inferred.
#[derive(Debug)]
pub struct MaterialSlabLaw {
    card: Arc<MaterialCard>,
    original: SpectralOpticalStatePoint,
    temperature_axis: String,
    temperature_bounds_k: (f64, f64),
    max_knot_span_nm: Option<f64>,
    identity: ContentHash,
}

impl MaterialSlabLaw {
    /// Retain an original dielectric source and its typed temperature axis.
    ///
    /// Tabulated sources require `Some` positive finite phase-index knot gap;
    /// analytic Sellmeier sources require `None`. The declared temperature
    /// interval further restricts every evaluation; it does not certify source
    /// coverage between its endpoints. Each accepted state must independently
    /// satisfy every pinned claim and the original optical owner's admission.
    /// `max_source_properties` bounds the complete original property bundle,
    /// which native consumers charge once per distinct law in their budgets.
    pub fn try_new(
        card: Arc<MaterialCard>,
        optical: SpectralOpticalStatePoint,
        temperature_axis: &str,
        temperature_bounds_k: (f64, f64),
        max_knot_span_nm: Option<f64>,
        max_source_properties: usize,
    ) -> Result<Self, SurfaceHeatingError> {
        if optical.resolved().properties().len() > max_source_properties {
            return Err(SurfaceHeatingError::Budget);
        }
        match (&optical, max_knot_span_nm) {
            (SpectralOpticalStatePoint::TabulatedDielectric(_), Some(gap))
                if gap.is_finite() && gap > 0.0 => {}
            (SpectralOpticalStatePoint::DielectricSellmeier(_), None) => {}
            _ => {
                return Err(invalid(
                    "material slab requires a dielectric source and its original gap convention",
                ));
            }
        }
        let (minimum_k, maximum_k) = temperature_bounds_k;
        if !minimum_k.is_finite()
            || !maximum_k.is_finite()
            || minimum_k <= 0.0
            || maximum_k < minimum_k
        {
            return Err(invalid(
                "material slab requires positive finite temperature support",
            ));
        }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let point = optical.resolved().query();
        if point.axis_quantities().get(temperature_axis) != Some(&absolute) {
            return Err(invalid(
                "material slab requires an existing typed absolute-temperature axis",
            ));
        }
        let initial = point
            .axes()
            .get(temperature_axis)
            .copied()
            .filter(|t| t.is_finite() && (minimum_k..=maximum_k).contains(t))
            .ok_or_else(|| invalid("original slab temperature must lie in the declared support"))?;
        let mut hash = DomainHasher::new("org.frankensim.render.material-temperature-slab.v1");
        hash.update(optical.identity().as_bytes());
        hash.update(&(temperature_axis.len() as u64).to_le_bytes());
        hash.update(temperature_axis.as_bytes());
        for value in [minimum_k, maximum_k] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.update(&[u8::from(max_knot_span_nm.is_some())]);
        if let Some(gap) = max_knot_span_nm {
            hash.update(&gap.to_bits().to_le_bytes());
        }
        let law = Self {
            card,
            original: optical,
            temperature_axis: temperature_axis.into(),
            temperature_bounds_k,
            max_knot_span_nm,
            identity: hash.finalize(),
        };
        // Validate original card/pins, all auxiliary claims, native optics and
        // spectral support before publishing the immutable temperature law.
        law.at_temperature(initial)?;
        Ok(law)
    }

    /// Complete original source, including its retained selected-claim receipts.
    pub const fn material_source(&self) -> &SpectralOpticalStatePoint {
        &self.original
    }

    pub const fn temperature_bounds_k(&self) -> (f64, f64) {
        self.temperature_bounds_k
    }

    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    pub(super) fn budget_size(&self) -> usize {
        self.original.resolved().properties().len()
    }

    fn source_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<SpectralOpticalStatePoint, SurfaceHeatingError> {
        let (lo, hi) = self.temperature_bounds_k;
        if !temperature_k.is_finite() || temperature_k < lo || temperature_k > hi {
            return Err(invalid(
                "material slab temperature lies outside its declared support",
            ));
        }
        let point = self
            .original
            .resolved()
            .query()
            .with(&self.temperature_axis, temperature_k)
            .map_err(|source| {
                SurfaceHeatingError::MaterialSource(MaterialStatePointError::Query {
                    property: format!("axis:{}", self.temperature_axis),
                    source,
                })
            })?;
        self.original
            .at_query_point(&self.card, &point)
            .map_err(SurfaceHeatingError::MaterialSource)
    }

    /// Resolve one smooth, bare optical snapshot at the requested temperature.
    /// It can feed the existing tracer or native slab emission; later spectral
    /// and directional queries reuse this snapshot without material queries.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<MaterialOpticalBinding, SurfaceHeatingError> {
        let source = self.source_at_temperature(temperature_k)?;
        let binding = match &source {
            SpectralOpticalStatePoint::TabulatedDielectric(optical) => {
                MaterialOpticalBinding::try_spectral_tabulated_dielectric(
                    optical.resolved(),
                    optical,
                    self.max_knot_span_nm
                        .ok_or_else(|| invalid("missing slab index gap"))?,
                    None,
                    self.identity,
                )
            }
            SpectralOpticalStatePoint::DielectricSellmeier(optical) => {
                MaterialOpticalBinding::try_spectral_sellmeier(
                    optical.resolved(),
                    optical,
                    None,
                    self.identity,
                )
            }
            SpectralOpticalStatePoint::Conductor(_) => {
                return Err(invalid(
                    "a material slab must retain its dielectric source family",
                ));
            }
        };
        binding.map_err(SurfaceHeatingError::MaterialOptics)
    }
}

#[cfg(test)]
mod tests;
