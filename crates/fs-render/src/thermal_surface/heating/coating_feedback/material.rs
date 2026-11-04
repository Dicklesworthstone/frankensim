//! Original finite-film n/k claims on the accepted native temperature clock.

use super::{Coordinate, Model, SurfaceCoatingLaw, SurfaceHeatingError, invalid};
use crate::conductor::{
    CONDUCTOR_IOR_SAMPLE_COUNT, CoatingIorSample, CoherentAbsorbingCoating, ConductorDataStatus,
    ConductorSource,
};
use fs_blake3::DomainHasher;
use fs_matdb::MaterialCard;
use fs_material::state_point::{MaterialStatePointError, SpectralConductorStatePoint};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use std::sync::Arc;

#[derive(Debug)]
pub(super) struct MaterialTemperatureCoating {
    card: Arc<MaterialCard>,
    original: SpectralConductorStatePoint,
    temperature_axis: String,
    thickness_m: f64,
    max_knot_span_nm: f64,
}

impl SurfaceCoatingLaw {
    /// Bind original n/k material claims to one uniform finite film.
    ///
    /// The original source name `SpectralConductorStatePoint` denotes its
    /// complex-index schema; that schema admits k=0. This finite-film consumer
    /// retains exact transparency and the original positive dissipative loss
    /// modes, without imposing the opaque substrate's k>0 requirement. The
    /// original coherent solver owns reflection, transmission, interference
    /// and absorption; only actual absorbed power heats the native receiver.
    ///
    /// Each accepted temperature refreshes every pinned claim, including
    /// source wavelengths and auxiliary properties. No temperature fit, RGB
    /// lift, alpha-to-k fit, changed claim selection or extrapolation is added.
    /// `temperature_axis` must already denote a typed static absolute
    /// temperature. All other coordinates and quantity conventions are held.
    /// The declared temperature interval is an additional restriction, not
    /// proof of continuous support by the source: every queried state may fail.
    ///
    /// Thickness and locally parallel geometry stay fixed. The film is
    /// isothermal in its assigned native cell, with prompt thermalization and
    /// coherent optics. No expansion, depth-resolved heat, phase mixture or
    /// thermal material/card association is inferred. Existing transmitting
    /// and opaque heating, emission and exchange paths consume this same law.
    ///
    /// `max_source_properties` bounds the complete original property bundle.
    /// Consumer optical-knot budgets charge that count once per distinct law.
    #[allow(clippy::too_many_arguments)]
    pub fn try_material_temperature(
        card: Arc<MaterialCard>,
        optical: SpectralConductorStatePoint,
        temperature_axis: &str,
        temperature_bounds_k: (f64, f64),
        thickness_m: f64,
        max_knot_span_nm: f64,
        max_source_properties: usize,
    ) -> Result<Self, SurfaceHeatingError> {
        if optical.resolved().properties().len() > max_source_properties {
            return Err(SurfaceHeatingError::Budget);
        }
        let (minimum_k, maximum_k) = temperature_bounds_k;
        if !minimum_k.is_finite()
            || !maximum_k.is_finite()
            || minimum_k <= 0.0
            || maximum_k < minimum_k
        {
            return Err(invalid(
                "material film requires positive finite temperature support",
            ));
        }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let point = optical.resolved().query();
        if point.axis_quantities().get(temperature_axis) != Some(&absolute) {
            return Err(invalid(
                "material film requires an existing typed absolute-temperature axis",
            ));
        }
        let initial = point
            .axes()
            .get(temperature_axis)
            .copied()
            .filter(|t| t.is_finite() && (minimum_k..=maximum_k).contains(t))
            .ok_or_else(|| invalid("original film temperature must lie in the declared support"))?;
        let mut hash = DomainHasher::new("org.frankensim.render.material-temperature-coating.v1");
        hash.update(optical.identity().as_bytes());
        hash.update(&(temperature_axis.len() as u64).to_le_bytes());
        hash.update(temperature_axis.as_bytes());
        for value in [minimum_k, maximum_k, thickness_m, max_knot_span_nm] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        let source = MaterialTemperatureCoating {
            card,
            original: optical,
            temperature_axis: temperature_axis.into(),
            thickness_m,
            max_knot_span_nm,
        };
        // Original card/pins, all auxiliary claims, optical support and native
        // thickness/index/gap restrictions are admitted before publication.
        source.coating_at_temperature(initial)?;
        Ok(Self {
            identity: hash.finalize(),
            coordinate: Coordinate::Temperature,
            temperature_bounds_k,
            model: Model::Material(Arc::new(source)),
        })
    }

    /// Original complete n/k source and its retained selected-claim receipts.
    /// Authored state rows have no original material bundle and return `None`.
    pub fn material_source(&self) -> Option<&SpectralConductorStatePoint> {
        match &self.model {
            Model::Material(source) => Some(&source.original),
            Model::Rows { .. } => None,
        }
    }
}

impl MaterialTemperatureCoating {
    pub(super) fn property_count(&self) -> usize {
        self.original.resolved().properties().len()
    }

    fn source_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<SpectralConductorStatePoint, SurfaceHeatingError> {
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

    pub(super) fn coating_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<CoherentAbsorbingCoating, SurfaceHeatingError> {
        let source = self.source_at_temperature(temperature_k)?;
        let original = source.samples();
        if original.len() != CONDUCTOR_IOR_SAMPLE_COUNT {
            return Err(invalid(
                "material film requires the original nine optical knots",
            ));
        }
        let first = original[0];
        let mut samples =
            [CoatingIorSample::try_new_in_band(first.wavelength_nm, first.eta, first.k)?;
                CONDUCTOR_IOR_SAMPLE_COUNT];
        for (sample, value) in samples.iter_mut().zip(original) {
            *sample = CoatingIorSample::try_new_in_band(value.wavelength_nm, value.eta, value.k)?;
        }
        Ok(CoherentAbsorbingCoating::try_new_in_band(
            samples,
            self.thickness_m,
            self.max_knot_span_nm,
            ConductorSource::try_new(source.identity(), ConductorDataStatus::MaterialStateClaims)?,
        )?)
    }
}

#[cfg(test)]
mod tests;
