//! Original dielectric sources on the accepted native slab temperature clock.

use super::{SurfaceHeatingError, invalid};
use crate::material_state::MaterialOpticalBinding;
use crate::material_state::thermo_optic::binding::MaterialThermoOpticBinding;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::state_point::{
    MaterialStatePointError, ResolvedMaterialStatePoint, SpectralDielectricAbsorptionStatePoint,
    SpectralOpticalStatePoint,
};
use fs_material::thermo_optic::SchottThermoOpticLaw;
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
/// Alternatively, an explicit SCHOTT mode retains phase and absorption sources
/// at the thermal law's reference temperature and refreshes only the separately
/// supplied mechanical state. Its absorption is deliberately fixed at reference.
/// An independently sourced absorption mode refreshes the original loss claims
/// at that same accepted state while keeping phase measurements at reference.
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
    thermo_optic: Option<MaterialThermoOpticBinding>,
    absorption_source: Option<SpectralDielectricAbsorptionStatePoint>,
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
            thermo_optic: None,
            absorption_source: None,
            identity: hash.finalize(),
        };
        // Validate original card/pins, all auxiliary claims, native optics and
        // spectral support before publishing the immutable temperature law.
        law.at_temperature(initial)?;
        Ok(law)
    }

    /// Use the shared SCHOTT phase law on the accepted native temperature clock,
    /// with explicitly frozen reference absorption and smooth bare surfaces.
    ///
    /// `optical` supplies the original phase/absorption claims at the law's T0.
    /// `reference_state` is the independently resolved physical state of that
    /// same card and complete condition. Its selected requirements are retained
    /// and re-resolved at every actual receiver temperature. A narrow mechanical
    /// bundle therefore permits optical measurements valid only at T0; no
    /// reference requirement is silently removed to extend source support.
    ///
    /// The original phase law is corrected exactly once. Later temperature
    /// dependence of the source's phase or absorption claims is not queried in
    /// this mode. The SCHOTT interval bounds thermal queries; mechanical source
    /// validity is independently enforced at each endpoint. No density-driven
    /// geometry change or phase-mixture model is inferred.
    ///
    /// Tabulated/Sellmeier gap conventions match [`Self::try_new`]. The property
    /// budget charges the sum of optical and mechanical reference requirements,
    /// conservatively retaining duplicates; native consumers charge a shared
    /// law once. Spectral/angular quadrature reuses the prepared native glass.
    pub fn try_schott_with_reference_absorption(
        card: Arc<MaterialCard>,
        optical: SpectralOpticalStatePoint,
        reference_state: &ResolvedMaterialStatePoint,
        law: SchottThermoOpticLaw,
        temperature_axis: &str,
        max_knot_span_nm: Option<f64>,
        max_source_properties: usize,
    ) -> Result<Self, SurfaceHeatingError> {
        let requirements = optical
            .resolved()
            .properties()
            .len()
            .checked_add(reference_state.properties().len())
            .ok_or(SurfaceHeatingError::Budget)?;
        if requirements > max_source_properties {
            return Err(SurfaceHeatingError::Budget);
        }
        let [minimum_k, maximum_k] = law.temperature_range_k();
        let mut source = Self::try_new(
            card,
            optical,
            temperature_axis,
            (minimum_k, maximum_k),
            max_knot_span_nm,
            max_source_properties,
        )?;
        let reference_binding = source.bind_source(reference_state, &source.original)?;
        let thermal = MaterialThermoOpticBinding::try_new(
            reference_binding,
            reference_state,
            law,
            temperature_axis,
        )
        .map_err(SurfaceHeatingError::MaterialThermoOptic)?;
        let mut hash = DomainHasher::new("org.frankensim.render.material-schott-reference-slab.v1");
        hash.update(source.identity.as_bytes());
        hash.update(reference_binding.identity().as_bytes());
        hash.update(law.identity().as_bytes());
        source.identity = hash.finalize();
        source.thermo_optic = Some(thermal);
        // Revalidate the selected physical requirements and their pinned claims
        // on this exact card, even when the original optics are only valid at T0.
        source.at_temperature(law.reference_temperature_k())?;
        Ok(source)
    }

    /// Use original SCHOTT refraction with separately sourced current absorption.
    ///
    /// Phase-index measurements remain at T0. The independent absorption source
    /// must start at that same complete reference condition; its selected loss,
    /// wavelength and auxiliary claims are then refreshed at each accepted
    /// receiver temperature. No missing phase claims, temperature fit or optical
    /// tail is introduced. Reference and current absorption are distinct modes.
    ///
    /// The property ceiling charges all three complete source bundles: original
    /// optics, actual mechanical requirements and independent absorption. All
    /// native illumination, emission and exchange consumers share this law and
    /// its ordinary prepared material binding without extra wavelength queries.
    #[allow(clippy::too_many_arguments)]
    pub fn try_schott_with_absorption_source(
        card: Arc<MaterialCard>,
        optical: SpectralOpticalStatePoint,
        reference_state: &ResolvedMaterialStatePoint,
        absorption: SpectralDielectricAbsorptionStatePoint,
        law: SchottThermoOpticLaw,
        temperature_axis: &str,
        max_knot_span_nm: Option<f64>,
        max_source_properties: usize,
    ) -> Result<Self, SurfaceHeatingError> {
        let requirements = optical
            .resolved()
            .properties()
            .len()
            .checked_add(reference_state.properties().len())
            .and_then(|count| count.checked_add(absorption.resolved().properties().len()))
            .ok_or(SurfaceHeatingError::Budget)?;
        if requirements > max_source_properties {
            return Err(SurfaceHeatingError::Budget);
        }
        let mut source = Self::try_schott_with_reference_absorption(
            card,
            optical,
            reference_state,
            law,
            temperature_axis,
            max_knot_span_nm,
            max_source_properties,
        )?;
        source
            .thermo_optic
            .as_ref()
            .expect("the SCHOTT constructor retains its original thermal binding")
            .at_material_state_with_absorption(reference_state, &absorption)
            .map_err(SurfaceHeatingError::MaterialThermoOptic)?;
        let mut hash = DomainHasher::new("org.frankensim.render.material-schott-absorbing-slab.v1");
        hash.update(source.identity.as_bytes());
        hash.update(absorption.identity().as_bytes());
        source.identity = hash.finalize();
        source.absorption_source = Some(absorption);
        // Re-resolve all original absorption pins against this exact card before
        // exposing the law, even when its initial material values already match.
        source.at_temperature(law.reference_temperature_k())?;
        Ok(source)
    }

    /// Complete original source, including its retained selected-claim receipts.
    pub const fn material_source(&self) -> &SpectralOpticalStatePoint {
        &self.original
    }

    /// Independent original loss source, when explicitly selected by the caller.
    pub const fn absorption_source(&self) -> Option<&SpectralDielectricAbsorptionStatePoint> {
        self.absorption_source.as_ref()
    }

    pub const fn temperature_bounds_k(&self) -> (f64, f64) {
        self.temperature_bounds_k
    }

    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    pub(super) fn budget_size(&self) -> usize {
        self.original.resolved().properties().len()
            + self
                .thermo_optic
                .as_ref()
                .map_or(0, |source| source.reference_state().properties().len())
            + self
                .absorption_source
                .as_ref()
                .map_or(0, |source| source.resolved().properties().len())
    }

    /// Read the actual native receiver without projecting or fitting its state.
    pub(super) fn coordinate_at(
        &self,
        thermal: &fs_thermal::enthalpy::DepositionSimulation<'_>,
        cell: fs_thermal::enthalpy::CellId,
        material: fs_thermal::enthalpy::MaterialId,
    ) -> Result<f64, SurfaceHeatingError> {
        let native = thermal
            .cells()
            .get(&cell)
            .ok_or_else(|| invalid("material slab receiver must be an active native cell"))?;
        if native.material_id != material {
            return Err(invalid("material slab receiver material changed"));
        }
        let temperature_k = thermal.temperature_k(cell)?;
        let (minimum_k, maximum_k) = self.temperature_bounds_k;
        if !temperature_k.is_finite() || !(minimum_k..=maximum_k).contains(&temperature_k) {
            return Err(SurfaceHeatingError::TemperatureDomain {
                thermal_cell: cell,
                temperature_k,
                minimum_k,
                maximum_k,
            });
        }
        Ok(temperature_k)
    }

    fn source_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<SpectralOpticalStatePoint, SurfaceHeatingError> {
        let point = self.query_at_temperature(self.original.resolved(), temperature_k)?;
        self.original
            .at_query_point(&self.card, &point)
            .map_err(SurfaceHeatingError::MaterialSource)
    }

    fn query_at_temperature(
        &self,
        original: &ResolvedMaterialStatePoint,
        temperature_k: f64,
    ) -> Result<QueryPoint, SurfaceHeatingError> {
        let (lo, hi) = self.temperature_bounds_k;
        if !temperature_k.is_finite() || temperature_k < lo || temperature_k > hi {
            return Err(invalid(
                "material slab temperature lies outside its declared support",
            ));
        }
        original
            .query()
            .with(&self.temperature_axis, temperature_k)
            .map_err(|source| {
                SurfaceHeatingError::MaterialSource(MaterialStatePointError::Query {
                    property: format!("axis:{}", self.temperature_axis),
                    source,
                })
            })
    }

    /// Resolve one smooth, bare optical snapshot at the requested temperature.
    /// It can feed the existing tracer or native slab emission; later spectral
    /// and directional queries reuse this snapshot without material queries.
    /// The explicitly selected SCHOTT mode refreshes mechanical requirements
    /// and corrects its original reference phase. Absorption stays at reference
    /// or refreshes its own original source, according to the selected mode.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<MaterialOpticalBinding, SurfaceHeatingError> {
        if let Some(source) = &self.thermo_optic {
            let reference = source.reference_state();
            let point = self.query_at_temperature(reference, temperature_k)?;
            let accepted = reference
                .at_query_point(&self.card, &point)
                .map_err(SurfaceHeatingError::MaterialSource)?;
            if let Some(absorption) = &self.absorption_source {
                let current = absorption
                    .at_query_point(&self.card, &point)
                    .map_err(SurfaceHeatingError::MaterialSource)?;
                return source
                    .at_material_state_with_absorption(&accepted, &current)
                    .map_err(SurfaceHeatingError::MaterialThermoOptic);
            }
            let observed = source
                .at_material_state(&accepted)
                .map_err(SurfaceHeatingError::MaterialThermoOptic)?;
            return Ok(*observed.binding());
        }
        let source = self.source_at_temperature(temperature_k)?;
        self.bind_source(source.resolved(), &source)
    }

    fn bind_source(
        &self,
        accepted: &ResolvedMaterialStatePoint,
        source: &SpectralOpticalStatePoint,
    ) -> Result<MaterialOpticalBinding, SurfaceHeatingError> {
        let binding = match source {
            SpectralOpticalStatePoint::TabulatedDielectric(optical) => {
                MaterialOpticalBinding::try_spectral_tabulated_dielectric(
                    accepted,
                    optical,
                    self.max_knot_span_nm
                        .ok_or_else(|| invalid("missing slab index gap"))?,
                    None,
                    self.identity,
                )
            }
            SpectralOpticalStatePoint::DielectricSellmeier(optical) => {
                MaterialOpticalBinding::try_spectral_sellmeier(
                    accepted,
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
