//! Original thermochemical primitives into analytic molecular optical jets.
//!
//! On one original linear-primitive span at common static pressure,
//! N_i = x_i N_A/V_m, N_i' = N (x_i' - x_i T'/T),
//! N_i'' = -2 (T'/T) N_i'. The original ideal-gas EOS supplies V_m;
//! temperature, active composition, model/phase admission and source receipts
//! remain those of AxialGasField::state_at. No calorically perfect surrogate,
//! finite-difference Hessian or independently authored density is introduced.

use super::*;
use crate::gradient_index::thermal::TemperatureField;
use fs_material::gas_field::{AxialGasField, AxialGasSampling};
use fs_material::gas_scattering::AVOGADRO_PER_MOL;
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

/// Admission or original-source refusal for a smooth thermochemical span.
#[derive(Clone, Debug, PartialEq)]
pub enum AxialGasCompositionError {
    UnsupportedReconstruction,
    UnsupportedPhase,
    InvalidSpan,
    InvalidEmbedding,
    MissingSnapshotIdentity,
    ModelMismatch,
    OutsideSpan,
    NumericalRange,
    Material(fs_material::MaterialError),
}
impl AxialGasCompositionError {
    fn into_gradient(self) -> GradientIndexError {
        match self {
            Self::OutsideSpan => GradientIndexError::OutsideFieldDomain,
            Self::NumericalRange => GradientIndexError::NumericalRange,
            _ => GradientIndexError::InvalidSample,
        }
    }
}
impl core::fmt::Display for AxialGasCompositionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "axial thermochemical optics: {self:?}")
    }
}
impl core::error::Error for AxialGasCompositionError {}

/// One borrowed smooth span of an original source-resolved gas field.
///
/// The optical catalogue must match every canonical original species and its
/// molar mass, including species absent at a span endpoint. Zero local species
/// remain zero; their analytic entrance derivatives are retained. Molecular
/// phase/group dispersion and elastic scattering use the existing optical
/// owners; true absorption remains exactly zero under that molecular law.
/// The same object implements [`TemperatureField`] from the original state.
/// Canonical species and exact molar masses establish the mapping; optical
/// references remain explicitly assigned independent sources. They do not
/// replace the field's original thermochemical/transport source authority.
/// Derivative coefficients and density bounds are prepared once. Each query
/// still calls the original `state_at`, including its allocations and source
/// validity checks; endpoint admission alone cannot establish the validity of
/// intermediate states. Wavelength preparation caches only optical spectra.
///
/// Linear primitive fields can have gradient jumps at station knots. This
/// binding admits ONLY one selected span, with one-sided endpoint derivatives,
/// and refuses any query outside it. It never supplies a fictitious smooth
/// Hessian across a knot. Physical ray paths and interface collars must stay
/// inside the selected source support; crossing other spans requires explicit
/// domain/geometry handling. World embedding declares radial uniformity along
/// the plane perpendicular to the supplied axial direction.
pub struct AxialGasCompositionField<'a> {
    field: &'a AxialGasField,
    optics: &'a MolecularRefractionModel,
    span: usize,
    bounds: [f64; 2],
    origin: Point3,
    axis: Vec3,
    temperature_slope: f64,
    fraction_slope_numerators: [f64; MAX_GAS_MIXTURE_COMPONENTS],
    upper: [f64; MAX_GAS_MIXTURE_COMPONENTS],
    identity: ContentHash,
}

impl<'a> AxialGasCompositionField<'a> {
    /// Bind a source span and its physical world embedding without resampling.
    /// `span_index` selects stations i and i+1. Only `LinearPrimitives` and
    /// `FrozenGas` are admitted; nonideal/equilibrating phases are not silently
    /// differentiated as an ideal gas. A finite nonzero direction is normalized.
    /// The nonzero identity declares the frozen original physical snapshot.
    pub fn try_new(
        field: &'a AxialGasField,
        optics: &'a MolecularRefractionModel,
        span_index: usize,
        origin_m: Point3,
        axial_direction: Vec3,
        snapshot_identity: ContentHash,
    ) -> Result<Self, AxialGasCompositionError> {
        use AxialGasCompositionError as Error;
        if field.sampling() != AxialGasSampling::LinearPrimitives {
            return Err(Error::UnsupportedReconstruction);
        }
        if field.phase() != ThermochemicalPhase::FrozenGas {
            return Err(Error::UnsupportedPhase);
        }
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(Error::MissingSnapshotIdentity);
        }
        if [
            origin_m.x,
            origin_m.y,
            origin_m.z,
            axial_direction.x,
            axial_direction.y,
            axial_direction.z,
        ]
        .into_iter()
        .any(|x| !x.is_finite())
        {
            return Err(Error::InvalidEmbedding);
        }
        let scale = axial_direction
            .x
            .abs()
            .max(axial_direction.y.abs())
            .max(axial_direction.z.abs());
        if scale == 0.0 {
            return Err(Error::InvalidEmbedding);
        }
        let scaled = Vec3::new(
            axial_direction.x / scale,
            axial_direction.y / scale,
            axial_direction.z / scale,
        );
        let axis = scaled.scale(1.0 / fs_math::det::sqrt(scaled.dot(scaled)));
        let end = span_index.checked_add(1).ok_or(Error::InvalidSpan)?;
        let bounds = [
            *field
                .stations_m()
                .get(span_index)
                .ok_or(Error::InvalidSpan)?,
            *field.stations_m().get(end).ok_or(Error::InvalidSpan)?,
        ];
        let original = field.model().thermodynamics().components();
        if optics.components().len() != original.len()
            || original.len() > MAX_GAS_MIXTURE_COMPONENTS
            || optics
                .components()
                .iter()
                .zip(original)
                .any(|(optical, species)| {
                    optical.species() != Some(species.species().as_str())
                        || optical.molar_mass_kg_mol() != species.molar_mass().value()
                })
        {
            return Err(Error::ModelMismatch);
        }
        let left = &field.primitive_states()[span_index];
        let right = &field.primitive_states()[end];
        let width = bounds[1] - bounds[0];
        let temperature_slope = (right.temperature_k - left.temperature_k) / width;
        let mut fraction_slope_numerators = [0.0; MAX_GAS_MIXTURE_COMPONENTS];
        let mut hash = DomainHasher::new("org.frankensim.render.axial-gas-composition.v1");
        hash.update(snapshot_identity.as_bytes());
        hash.update(optics.content_identity().as_bytes());
        hash.update(&(span_index as u64).to_le_bytes());
        for value in [
            origin_m.x,
            origin_m.y,
            origin_m.z,
            axis.x,
            axis.y,
            axis.z,
            bounds[0],
            bounds[1],
            field.static_pressure_pa(),
            left.temperature_k,
            right.temperature_k,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        for (i, (&a, &b)) in left
            .composition
            .fractions()
            .iter()
            .zip(right.composition.fractions())
            .enumerate()
        {
            // The rational derivative's numerator is constant throughout the
            // span. Prepare it once, avoiding cancellation against rounded
            // interpolated fractions at every optical stage.
            fraction_slope_numerators[i] =
                ((b - a) / width) * left.temperature_k - a * temperature_slope;
            hash.update(&a.to_bits().to_le_bytes());
            hash.update(&b.to_bits().to_le_bytes());
        }
        if !temperature_slope.is_finite()
            || fraction_slope_numerators.iter().any(|x| !x.is_finite())
        {
            return Err(Error::NumericalRange);
        }
        let mut result = Self {
            field,
            optics,
            span: span_index,
            bounds,
            origin: origin_m,
            axis,
            temperature_slope,
            fraction_slope_numerators,
            upper: [0.0; MAX_GAS_MIXTURE_COMPONENTS],
            identity: hash.finalize(),
        };
        // A positive-temperature linear/linear ratio is monotone on this span.
        // Preserve original endpoint values, with outward rounding. This is a
        // numerical envelope, not an interval certificate; tracking still
        // checks every sampled value against the retained bound.
        let mut upper = [0.0_f64; MAX_GAS_MIXTURE_COMPONENTS];
        for x in bounds {
            let state = field.state_at(x).map_err(Error::Material)?;
            result.with_jets(&state, |i, density, _, _| {
                upper[i] = upper[i].max(density);
            })?;
        }
        for value in &mut upper[..original.len()] {
            if *value > 0.0 {
                *value = value.next_up();
            }
            if !value.is_finite() {
                return Err(Error::NumericalRange);
            }
        }
        result.upper = upper;
        Ok(result)
    }

    /// The retained thermochemical field and its unmodified source models.
    pub const fn source(&self) -> &'a AxialGasField {
        self.field
    }
    /// Index of the original pair of stations defining this smooth domain.
    pub const fn span_index(&self) -> usize {
        self.span
    }
    /// Closed axial support in metres, measured from the supplied world origin.
    pub const fn axial_bounds_m(&self) -> [f64; 2] {
        self.bounds
    }
    /// Physical snapshot, original primitives, optical source and embedding.
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }
    /// Source-derived per-species number-density envelope in molecules/m³.
    /// Use directly with the existing `bind_scattering_field`; its original
    /// spectral/dilute admission and runtime bound verification still apply.
    pub fn component_number_density_upper_bounds(&self) -> &[f64] {
        &self.upper[..self.optics.components().len()]
    }
    /// Feed this exact source into ordinary, wavelength-prepared and second-
    /// order native refractive transport through the existing molecular owner.
    pub fn bind_refraction(
        &self,
    ) -> Result<MaterialGasCompositionRefractiveField<'_, Self>, GasRefractiveFieldError> {
        MaterialGasCompositionRefractiveField::try_new(self.optics, self, self.identity)
    }
    /// Original accepted local gas state, including source/phase diagnostics.
    /// No phase, EOS, optical spectrum or thermal state is reconstructed here.
    pub fn sample_state(
        &self,
        position_m: Point3,
    ) -> Result<ResolvedThermochemicalGasState, AxialGasCompositionError> {
        let x = position_m.delta_from(self.origin).dot(self.axis);
        if !x.is_finite() || x < self.bounds[0] || x > self.bounds[1] {
            return Err(AxialGasCompositionError::OutsideSpan);
        }
        self.field
            .state_at(x)
            .map_err(AxialGasCompositionError::Material)
    }

    fn with_jets(
        &self,
        state: &ResolvedThermochemicalGasState,
        mut emit: impl FnMut(usize, f64, f64, f64),
    ) -> Result<(), AxialGasCompositionError> {
        let number = AVOGADRO_PER_MOL / state.equation_of_state().molar_volume().value();
        let ratio = self.temperature_slope / state.state().temperature;
        if !number.is_finite() || number <= 0.0 || !ratio.is_finite() {
            return Err(AxialGasCompositionError::NumericalRange);
        }
        let receipt = state.thermodynamics().receipt();
        let mut active = 0;
        for (i, component) in self.optics.components().iter().enumerate() {
            let fraction = if receipt
                .species()
                .get(active)
                .is_some_and(|name| component.species() == Some(name.as_str()))
            {
                let fraction = f64::from_bits(receipt.mole_fraction_bits()[active]);
                active += 1;
                fraction
            } else {
                0.0
            };
            let density = number * fraction;
            let slope = number * (self.fraction_slope_numerators[i] / state.state().temperature);
            let curvature = -2.0 * (ratio * slope);
            if [density, slope, curvature]
                .into_iter()
                .any(|x| !x.is_finite())
                || (fraction > 0.0 && density == 0.0)
            {
                return Err(AxialGasCompositionError::NumericalRange);
            }
            emit(i, density, slope, curvature);
        }
        if active != receipt.species().len() {
            return Err(AxialGasCompositionError::ModelMismatch);
        }
        Ok(())
    }
}

impl GasCompositionGradientField for AxialGasCompositionField<'_> {
    fn model_identity(&self) -> ContentHash {
        self.optics.content_identity()
    }
    fn component_count(&self) -> usize {
        self.optics.components().len()
    }
    fn density_basis(&self) -> PartialDensityBasis {
        PartialDensityBasis::MoleculesPerM3
    }
    fn sample_partial_densities(
        &self,
        p: Point3,
        output: &mut [GasPartialDensitySample],
    ) -> Result<(), GradientIndexError> {
        if output.len() != self.component_count() {
            return Err(GradientIndexError::InvalidSample);
        }
        let state = self
            .sample_state(p)
            .map_err(AxialGasCompositionError::into_gradient)?;
        self.with_jets(&state, |i, density, slope, _| {
            output[i] = GasPartialDensitySample {
                density,
                gradient_per_m: self.axis.scale(slope),
            };
        })
        .map_err(AxialGasCompositionError::into_gradient)
    }
}
impl GasCompositionSecondOrderField for AxialGasCompositionField<'_> {
    fn sample_partial_density_jets(
        &self,
        p: Point3,
        output: &mut [GasPartialDensityJet],
    ) -> Result<(), GradientIndexError> {
        if output.len() != self.component_count() {
            return Err(GradientIndexError::InvalidSample);
        }
        let state = self
            .sample_state(p)
            .map_err(AxialGasCompositionError::into_gradient)?;
        let [x, y, z] = [self.axis.x, self.axis.y, self.axis.z];
        self.with_jets(&state, |i, density, slope, curvature| {
            output[i] = GasPartialDensityJet {
                sample: GasPartialDensitySample {
                    density,
                    gradient_per_m: self.axis.scale(slope),
                },
                hessian_diagonal_per_m2: [curvature * x * x, curvature * y * y, curvature * z * z],
                hessian_off_diagonal_per_m2: [
                    curvature * x * y,
                    curvature * x * z,
                    curvature * y * z,
                ],
            };
        })
        .map_err(AxialGasCompositionError::into_gradient)
    }
}
impl TemperatureField for AxialGasCompositionField<'_> {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        self.sample_state(p)
            .map(|state| state.state().temperature)
            .map_err(AxialGasCompositionError::into_gradient)
    }
}

#[cfg(test)]
mod tests;
