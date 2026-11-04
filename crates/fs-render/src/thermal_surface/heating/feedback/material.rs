//! Original material-card claims evaluated on the native temperature clock.

use super::{Coordinate, Model, SurfaceConductorLaw, SurfaceHeatingError, invalid};
use crate::conductor::ConductorOptics;
use crate::material_state::{MaterialOpticalBinding, OpticalAppearance};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::MaterialCard;
use fs_material::state_point::{MaterialStatePointError, SpectralConductorStatePoint};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use std::sync::Arc;

#[derive(Debug)]
pub(super) struct MaterialTemperatureConductor {
    card: Arc<MaterialCard>,
    original: SpectralConductorStatePoint,
    temperature_axis: String,
    max_knot_span_nm: f64,
    identity: ContentHash,
}

impl SurfaceConductorLaw {
    /// Follow the original pinned optical claims at one native temperature.
    ///
    /// `temperature_axis` must already be an explicitly typed static absolute
    /// temperature in the complete original bundle. Every other coordinate and
    /// quantity convention stays fixed. Wavelengths, n/k, density and auxiliary
    /// properties are refreshed together through the original selected claims;
    /// no new temperature fit or replacement-claim selection is introduced.
    ///
    /// The temperature bounds are additional caller restrictions, not a claim
    /// that the material supports every interior point. Each native endpoint
    /// is queried and admitted, including dark endpoints. Unsupported sources
    /// reject the thermal candidate before publication. This explicit feedback
    /// still requires timestep refinement; it is not an implicit optical solve.
    ///
    /// This is the original smooth, opaque, isothermal, semi-infinite substrate
    /// model. Native conductor admission requires strictly positive n and k at
    /// every source knot; exact transparency is refused rather than assigned
    /// heat. A finite weakly absorbing specimen needs the existing transmitting
    /// volume/film model and its thickness. Neither deposition depth nor a
    /// thermal material/card association is inferred from the optical source.
    ///
    /// `max_source_properties` bounds all retained scalar queries, including
    /// auxiliary claims. Existing consumer optical-knot budgets charge that
    /// same property count once per distinct law. Source claims keep their
    /// original authority and finite wavelength support, including supplied IR.
    #[allow(clippy::too_many_arguments)]
    pub fn try_material_temperature(
        card: Arc<MaterialCard>,
        optical: SpectralConductorStatePoint,
        temperature_axis: &str,
        temperature_bounds_k: (f64, f64),
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
                "material optics require positive finite temperature support",
            ));
        }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let point = optical.resolved().query();
        if point.axis_quantities().get(temperature_axis) != Some(&absolute) {
            return Err(invalid(
                "material optical feedback requires an existing typed absolute-temperature axis",
            ));
        }
        let initial = point
            .axes()
            .get(temperature_axis)
            .copied()
            .filter(|t| t.is_finite() && (minimum_k..=maximum_k).contains(t))
            .ok_or_else(|| {
                invalid("original optical temperature must lie in the declared support")
            })?;
        let mut hash = DomainHasher::new("org.frankensim.render.material-temperature-conductor.v1");
        hash.update(optical.identity().as_bytes());
        hash.update(&(temperature_axis.len() as u64).to_le_bytes());
        hash.update(temperature_axis.as_bytes());
        for value in [minimum_k, maximum_k, max_knot_span_nm] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        let identity = hash.finalize();
        let source = MaterialTemperatureConductor {
            card,
            original: optical,
            temperature_axis: temperature_axis.into(),
            max_knot_span_nm,
            identity,
        };
        // Validate the exact card, original pins, current source band, knot-gap
        // limit and opaque n/k domain before admitting any feedback law.
        source.optics_at_temperature(initial)?;
        Ok(Self {
            identity,
            coordinate: Coordinate::Temperature,
            temperature_bounds_k,
            model: Model::Material(Arc::new(source)),
        })
    }

    /// Original complete source bundle and pinned receipts, without any fitted
    /// temperature rows. Authored row laws return `None`.
    pub fn material_source(&self) -> Option<&SpectralConductorStatePoint> {
        match &self.model {
            Model::Material(source) => Some(&source.original),
            Model::Rows { .. } => None,
        }
    }
}

impl MaterialTemperatureConductor {
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

    pub(super) fn optics_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ConductorOptics, SurfaceHeatingError> {
        let source = self.source_at_temperature(temperature_k)?;
        // Reuse the material/surface ingress, including strict native opacity,
        // original source knots, full-state identity and finite-band admission.
        let binding = MaterialOpticalBinding::try_spectral_conductor(
            source.resolved(),
            &source,
            self.max_knot_span_nm,
            0.0,
            self.identity,
        )
        .map_err(SurfaceHeatingError::MaterialOptics)?;
        match binding.appearance() {
            OpticalAppearance::Conductor { optics, .. } => Ok(optics),
            OpticalAppearance::Dielectric { .. } => Err(invalid(
                "material conductor binding must remain a conductor",
            )),
        }
    }
}

#[cfg(test)]
mod tests;
