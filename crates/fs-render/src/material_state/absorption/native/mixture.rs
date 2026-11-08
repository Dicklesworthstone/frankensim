//! Homogeneous independent absorbers retain their original material sources.

use super::{
    AbsorbingGrid, CellPalette, ContentHash, Cx, DomainHasher, MaterialAbsorptionBinding,
    MaterialAbsorptionEmission, MaterialAbsorptionField, MaterialAbsorptionFieldError,
    MaterialAbsorptionPreparationBudget, PreparedMaterial, VolumeTransportError, checkpoint,
    planck_radiance_w_m2_sr_nm,
};
use crate::volumes::transport::{
    MAX_OPTICAL_MIXTURE_COMPONENTS, OpticalMixtureBudget, PositiveSum, weighted_product,
};
use fs_math::det;

/// One original accepted absorber and its explicit reference-mass loading.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialAbsorptionComponent<'a> {
    /// Source law, intrinsic density receipt and complete accepted condition.
    pub material: MaterialAbsorptionBinding<'a>,
    /// Kilograms of this constituent per kg of the rendered density's reference
    /// mass. Total-mixture fractions should sum to one; carrier loadings need
    /// not. This owner never normalizes them or infers a density basis.
    pub mass_weight: f64,
}

/// Original constituent bindings in one homogeneous independent-absorption law.
///
/// This borrows the declared component order without allocating or replacing
/// any source spectrum. Absorption adds as `kappa = sum(w_i*kappa_i)`; accepted
/// LTE emission adds as `rho*sum(w_i*kappa_i*B_lambda(T_i))`. Different component
/// temperatures remain distinct. This assumes additive bulk absorption and
/// does not describe interfaces, scattering, a homogenized refractive index,
/// subpixel coverage, thermal equilibration or chemical interactions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialAbsorptionMixture<'a> {
    components: &'a [MaterialAbsorptionComponent<'a>],
    source_bounds_nm: [f64; 2],
    identity: ContentHash,
}

impl<'a> MaterialAbsorptionMixture<'a> {
    /// Admit bounded, nonnegative mass weights and the original common source
    /// support. Zero-weight entries retain validation, support and identity;
    /// they create no absorption and do not require a temperature. An all-zero
    /// recipe is transparent. A nonempty recipe is required.
    pub fn try_new(
        cx: &Cx<'_>,
        components: &'a [MaterialAbsorptionComponent<'a>],
        budget: OpticalMixtureBudget,
    ) -> Result<Self, MaterialAbsorptionFieldError> {
        checkpoint(cx)?;
        if components.len() > budget.max_components
            || components.len() > MAX_OPTICAL_MIXTURE_COMPONENTS
        {
            return Err(VolumeTransportError::Budget.into());
        }
        if components.is_empty() {
            return Err(VolumeTransportError::Mixture.into());
        }
        let mut hash = DomainHasher::new("org.frankensim.render.material-absorption-mixture.v1");
        hash.update(&(components.len() as u64).to_le_bytes());
        let mut bounds = [0.0_f64, f64::INFINITY];
        for component in components {
            checkpoint(cx)?;
            if !component.mass_weight.is_finite() || component.mass_weight < 0.0 {
                return Err(VolumeTransportError::Mixture.into());
            }
            let support = component.material.wavelength_bounds_nm();
            bounds[0] = bounds[0].max(support[0]);
            bounds[1] = bounds[1].min(support[1]);
            hash.update(component.material.identity().as_bytes());
            hash.update(&component.mass_weight.to_bits().to_le_bytes());
        }
        if bounds[0] > bounds[1] {
            return Err(VolumeTransportError::Wavelength.into());
        }
        checkpoint(cx)?;
        Ok(Self {
            components,
            source_bounds_nm: bounds,
            identity: hash.finalize(),
        })
    }

    /// Original weights and accepted material bindings in declaration order.
    #[must_use]
    pub const fn components(&self) -> &'a [MaterialAbsorptionComponent<'a>] {
        self.components
    }

    /// Closed intersection of every declared source support, including zero
    /// weights. A prepared field evaluates only its immutable wavelength nodes.
    #[must_use]
    pub const fn source_wavelength_bounds_nm(&self) -> [f64; 2] {
        self.source_bounds_nm
    }

    /// Original binding identities, declaration order and exact mass weights.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    fn prepare<const N: usize>(
        &self,
        cx: &Cx<'_>,
        wavelengths_nm: &[f64; N],
        emission: MaterialAbsorptionEmission,
        model: usize,
    ) -> Result<PreparedMaterial<N>, MaterialAbsorptionFieldError> {
        // Bounded temporary storage; final storage remains one packet/recipe.
        // No constituent loop, material query or Planck evaluation reaches rays.
        let mut weighted = [[0.0; N]; MAX_OPTICAL_MIXTURE_COMPONENTS];
        let mut temperatures = [None; MAX_OPTICAL_MIXTURE_COMPONENTS];
        let mut totals: [PositiveSum; N] = core::array::from_fn(|_| PositiveSum::default());
        let mut first_temperature = None;
        let mut common_temperature = true;
        for (index, component) in self.components.iter().enumerate() {
            checkpoint(cx)?;
            let samples = component.material.sample(wavelengths_nm)?;
            for lane in 0..N {
                weighted[index][lane] =
                    weighted_product(component.mass_weight, samples.mass_absorption_m2_kg()[lane])?;
                totals[lane].add(weighted[index][lane])?;
            }
            if emission == MaterialAbsorptionEmission::AcceptedLte
                && weighted[index].iter().any(|value| *value > 0.0)
            {
                let temperature = component.material.accepted_temperature_k().ok_or(
                    MaterialAbsorptionFieldError::ComponentTemperature {
                        model,
                        component: index,
                    },
                )?;
                temperatures[index] = Some(temperature);
                if let Some(first) = first_temperature {
                    common_temperature &= first == temperature;
                } else {
                    first_temperature = Some(temperature);
                }
            }
        }
        let mut absorption = [0.0; N];
        for (lane, total) in totals.into_iter().enumerate() {
            absorption[lane] = total.finish()?;
        }
        let source = if let Some(first) = first_temperature {
            let mut source = [0.0; N];
            for lane in 0..N {
                checkpoint(cx)?;
                if absorption[lane] == 0.0 {
                    continue;
                }
                if common_temperature {
                    source[lane] = planck_radiance_w_m2_sr_nm(first, wavelengths_nm[lane])
                        .map_err(|_| VolumeTransportError::Unrepresentable)?;
                    continue;
                }
                let mut radiance = PositiveSum::default();
                for index in 0..self.components.len() {
                    let coefficient = weighted[index][lane];
                    if coefficient == 0.0 {
                        continue;
                    }
                    let component_source = planck_radiance_w_m2_sr_nm(
                        temperatures[index].expect("admitted emitting constituent"),
                        wavelengths_nm[lane],
                    )
                    .map_err(|_| VolumeTransportError::Unrepresentable)?;
                    // Divide before multiplying so kappa*B need not fit. Match
                    // the original constituent thermal owner's range fallback.
                    let mut contribution = (coefficient / absorption[lane]) * component_source;
                    if contribution == 0.0 && component_source > 0.0 {
                        contribution = det::exp(
                            det::ln(coefficient) - det::ln(absorption[lane])
                                + det::ln(component_source),
                        );
                    }
                    radiance.add(contribution)?;
                }
                source[lane] = radiance.finish()?;
            }
            Some(source)
        } else {
            None
        };
        Ok(PreparedMaterial { absorption, source })
    }
}

/// Constituent recipes using the same prepared field, cameras and beam-energy
/// transport as individual materials. Only preparation knows the constituents.
pub type MaterialAbsorptionMixtureField<'g, 'field, 'material, const N: usize> =
    MaterialAbsorptionField<'g, 'field, 'material, N, MaterialAbsorptionMixture<'material>>;

impl<'g, 'field, 'material, const N: usize>
    MaterialAbsorptionField<'g, 'field, 'material, N, MaterialAbsorptionMixture<'material>>
{
    /// Evaluate original constituent laws at the exact nodes and prepare each
    /// homogeneous recipe once. The grid density is the explicitly declared
    /// reference mass per bulk m³. Component count times node count is charged
    /// for every recipe before allocation, including unused and zero weights.
    /// Output has no off-node interpolation or per-ray constituent work.
    pub fn try_new_mixtures(
        cx: &Cx<'_>,
        geometry: &'g AbsorbingGrid<'field>,
        materials: CellPalette<MaterialAbsorptionMixture<'material>>,
        wavelengths_nm: [f64; N],
        emission: MaterialAbsorptionEmission,
        budget: MaterialAbsorptionPreparationBudget,
    ) -> Result<Self, MaterialAbsorptionFieldError> {
        Self::prepare(
            cx,
            geometry,
            materials,
            wavelengths_nm,
            emission,
            budget,
            "org.frankensim.render.material-absorption-mixture-field.v1",
            |mixture| mixture.components.len(),
            |mixture, index| {
                Ok((
                    mixture.prepare(cx, &wavelengths_nm, emission, index)?,
                    mixture.source_bounds_nm,
                    mixture.identity,
                ))
            },
        )
    }
}

#[cfg(test)]
mod tests;
