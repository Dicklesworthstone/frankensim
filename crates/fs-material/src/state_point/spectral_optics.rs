//! Source-selected optical tables on their own finite vacuum-wavelength band.
//!
//! Wavelength coordinates are material properties in SI metres, so the grid,
//! coefficients, source claims and complete physical state share one identity.
//! These explicitly selected schemas do not change the legacy visible resolver.

use super::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, SpectralComplexIndexKnot,
    SpectralComplexIndexSample, SpectralComplexIndexStatePoint,
    VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY, resolve_material_state_point,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

/// Number of original knots retained by each finite-band optical table.
pub const SPECTRAL_OPTICAL_SAMPLE_COUNT: usize = 9;

/// Strictly increasing vacuum wavelengths in SI metres, in sample order.
/// These are source coordinates, not distances travelled inside the material.
pub const SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES: [&str; SPECTRAL_OPTICAL_SAMPLE_COUNT] = [
    "optical_spectral_wavelength_0_m",
    "optical_spectral_wavelength_1_m",
    "optical_spectral_wavelength_2_m",
    "optical_spectral_wavelength_3_m",
    "optical_spectral_wavelength_4_m",
    "optical_spectral_wavelength_5_m",
    "optical_spectral_wavelength_6_m",
    "optical_spectral_wavelength_7_m",
    "optical_spectral_wavelength_8_m",
];
/// Absolute real complex-index components at the corresponding source wavelengths.
pub const SPECTRAL_CONDUCTOR_ETA_PROPERTIES: [&str; SPECTRAL_OPTICAL_SAMPLE_COUNT] = [
    "optical_spectral_conductor_eta_0",
    "optical_spectral_conductor_eta_1",
    "optical_spectral_conductor_eta_2",
    "optical_spectral_conductor_eta_3",
    "optical_spectral_conductor_eta_4",
    "optical_spectral_conductor_eta_5",
    "optical_spectral_conductor_eta_6",
    "optical_spectral_conductor_eta_7",
    "optical_spectral_conductor_eta_8",
];
/// Nonnegative dimensionless extinction indices, not intensity absorption.
pub const SPECTRAL_CONDUCTOR_K_PROPERTIES: [&str; SPECTRAL_OPTICAL_SAMPLE_COUNT] = [
    "optical_spectral_conductor_k_0",
    "optical_spectral_conductor_k_1",
    "optical_spectral_conductor_k_2",
    "optical_spectral_conductor_k_3",
    "optical_spectral_conductor_k_4",
    "optical_spectral_conductor_k_5",
    "optical_spectral_conductor_k_6",
    "optical_spectral_conductor_k_7",
    "optical_spectral_conductor_k_8",
];
/// Dimensionless optical extinction index, distinct from attenuation per metre.
pub const SPECTRAL_CONDUCTOR_K_QUANTITY: QuantitySpec = QuantitySpec::semantic(SemanticType::new(
    QuantityKind::OpticalExtinctionCoefficient,
    ValueForm::Static,
));
/// Absolute real phase indices for the weakly absorbing dielectric model.
pub const SPECTRAL_DIELECTRIC_IOR_PROPERTIES: [&str; SPECTRAL_OPTICAL_SAMPLE_COUNT] = [
    "optical_spectral_dielectric_eta_0",
    "optical_spectral_dielectric_eta_1",
    "optical_spectral_dielectric_eta_2",
    "optical_spectral_dielectric_eta_3",
    "optical_spectral_dielectric_eta_4",
    "optical_spectral_dielectric_eta_5",
    "optical_spectral_dielectric_eta_6",
    "optical_spectral_dielectric_eta_7",
    "optical_spectral_dielectric_eta_8",
];
/// Internal intensity absorption at the source wavelengths [power Np/m].
/// Interface reflection and scattering losses are not part of these values.
pub const SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES: [&str; SPECTRAL_OPTICAL_SAMPLE_COUNT] = [
    "optical_spectral_absorption_per_m_0",
    "optical_spectral_absorption_per_m_1",
    "optical_spectral_absorption_per_m_2",
    "optical_spectral_absorption_per_m_3",
    "optical_spectral_absorption_per_m_4",
    "optical_spectral_absorption_per_m_5",
    "optical_spectral_absorption_per_m_6",
    "optical_spectral_absorption_per_m_7",
    "optical_spectral_absorption_per_m_8",
];
/// Exact quantity convention in `I(d) = I(0) * exp(-alpha * d)`.
pub const SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY: QuantitySpec =
    VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY;

const WAVELENGTH_QUANTITY: QuantitySpec = QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0]));
const INDEX_QUANTITY: QuantitySpec = QuantitySpec::dimensional(Dims::NONE);

/// One original real phase-index and intensity-absorption sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralDielectricSample {
    /// Vacuum wavelength [m].
    pub wavelength_m: f64,
    /// Absolute phase index in the native weakly absorbing domain `[1, 3]`.
    pub phase_index: f64,
    /// Internal intensity absorption [power Np/m].
    pub absorption_per_m: f64,
}

/// Source-wavelength ingress into the original complex-index material owner.
///
/// This wrapper resolves nine SI wavelength properties and their canonical n/k
/// roles, then delegates optical admission, interpolation and identity to
/// [`SpectralComplexIndexStatePoint`]. It adds no second complex-index model.
/// The same projection can drive surfaces, particle scattering and EM consumers.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralConductorStatePoint {
    complex_index: SpectralComplexIndexStatePoint,
    wavelengths_m: [f64; SPECTRAL_OPTICAL_SAMPLE_COUNT],
}

impl SpectralConductorStatePoint {
    /// Adopt a complete original snapshot without discarding any other property.
    /// All canonical optical quantities and values are rechecked, even when the
    /// original consumer used more permissive scalar admissibility requirements.
    pub fn try_from_material(
        resolved: ResolvedMaterialStatePoint,
    ) -> Result<Self, MaterialStatePointError> {
        let wavelengths_m = source_wavelengths(&resolved)?;
        let bindings: [SpectralComplexIndexKnot; SPECTRAL_OPTICAL_SAMPLE_COUNT] =
            core::array::from_fn(|index| SpectralComplexIndexKnot {
                wavelength_nm: wavelengths_m[index] * 1.0e9,
                real_index_property: PropertyKey::with_quantity(
                    SPECTRAL_CONDUCTOR_ETA_PROPERTIES[index],
                    INDEX_QUANTITY,
                ),
                extinction_index_property: PropertyKey::with_quantity(
                    SPECTRAL_CONDUCTOR_K_PROPERTIES[index],
                    SPECTRAL_CONDUCTOR_K_QUANTITY,
                ),
            });
        let complex_index = SpectralComplexIndexStatePoint::try_from_resolved(resolved, &bindings)?;
        Ok(Self {
            complex_index,
            wavelengths_m,
        })
    }

    /// Original card, complete state and every selected property-use receipt.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        self.complex_index.resolved()
    }

    /// The shared original n/k projection, including its source bindings and
    /// allocation-free spectral evaluation. No fitting or resampling is applied.
    /// Continue through this wrapper's `at_query_point` when SI wavelength
    /// properties evolve; a bare projection retains its authored bindings.
    #[must_use]
    pub const fn complex_index(&self) -> &SpectralComplexIndexStatePoint {
        &self.complex_index
    }

    /// The original complex-index owner's identity, including wavelength/key
    /// roles and the complete material bundle with SI wavelength properties.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.complex_index.identity()
    }

    /// The original owner's knots in vacuum nanometres, converted once from
    /// the selected SI wavelength properties. This ingress always has nine knots.
    #[must_use]
    pub fn samples(&self) -> &[SpectralComplexIndexSample] {
        self.complex_index.samples()
    }

    /// Closed source support in vacuum metres at this material state.
    #[must_use]
    pub fn wavelength_bounds_m(&self) -> (f64, f64) {
        (self.wavelengths_m[0], self.wavelengths_m[8])
    }

    /// Refresh the exact original claims and coordinate meanings atomically.
    /// Every refreshed wavelength and coefficient is readmitted; unsupported
    /// source states do not switch claims or mutate this snapshot.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        // Rebuild the bindings as well as values: these wavelengths are source
        // properties that may change with the accepted physical condition.
        Self::try_from_material(self.resolved().at_query_point(card, point)?)
    }
}

/// Real-index geometric optics and bulk absorption on one original source band.
///
/// This is an isotropic, weakly absorbing dielectric, not a complex-interface
/// model for strongly absorbing media. The samples alone establish neither
/// birefringence nor a broadband thermal-emissivity law.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralTabulatedDielectricStatePoint {
    resolved: ResolvedMaterialStatePoint,
    samples: [SpectralDielectricSample; SPECTRAL_OPTICAL_SAMPLE_COUNT],
}

impl SpectralTabulatedDielectricStatePoint {
    /// Adopt all original wavelength, index and typed absorption properties.
    /// Missing, ambiguous, noncanonical or nonphysical values refuse together.
    pub fn try_from_material(
        resolved: ResolvedMaterialStatePoint,
    ) -> Result<Self, MaterialStatePointError> {
        let wavelengths = source_wavelengths(&resolved)?;
        let mut samples = [SpectralDielectricSample {
            wavelength_m: 0.0,
            phase_index: 0.0,
            absorption_per_m: 0.0,
        }; 9];
        for (index, sample) in samples.iter_mut().enumerate() {
            let phase_index = source_value(
                &resolved,
                SPECTRAL_DIELECTRIC_IOR_PROPERTIES[index],
                INDEX_QUANTITY,
                ScalarAdmissibility::StrictlyPositive,
            )?;
            if !(1.0..=3.0).contains(&phase_index) {
                return Err(MaterialStatePointError::InvalidDerived {
                    quantity: "spectral_dielectric_phase_index",
                });
            }
            *sample = SpectralDielectricSample {
                wavelength_m: wavelengths[index],
                phase_index,
                absorption_per_m: source_value(
                    &resolved,
                    SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[index],
                    SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
                    ScalarAdmissibility::NonNegative,
                )?,
            };
        }
        Ok(Self { resolved, samples })
    }

    /// Original card, complete state and every selected property-use receipt.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Identity of this dielectric family and its complete original snapshot.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("org.frankensim.material.spectral-dielectric.v1");
        hash.update(self.resolved.identity().as_bytes());
        hash.finalize()
    }

    /// Original phase and absorption knots on their shared vacuum-wavelength grid.
    #[must_use]
    pub const fn samples(&self) -> &[SpectralDielectricSample; SPECTRAL_OPTICAL_SAMPLE_COUNT] {
        &self.samples
    }

    /// Closed source support in vacuum metres at this material state.
    #[must_use]
    pub fn wavelength_bounds_m(&self) -> (f64, f64) {
        (self.samples[0].wavelength_m, self.samples[8].wavelength_m)
    }

    /// Continue the original source claims at one new complete physical state.
    /// This does not certify the path between endpoints or evolve the specimen.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_material(self.resolved.at_query_point(card, point)?)
    }
}

/// Explicitly selected finite-band optical family. No card-name dispatch or
/// automatic competition with the legacy visible schemas is performed.
#[derive(Clone, Debug, PartialEq)]
pub enum SpectralOpticalStatePoint {
    /// Opaque complex-index conductor response.
    Conductor(SpectralConductorStatePoint),
    /// Weakly absorbing real-index dielectric response.
    TabulatedDielectric(SpectralTabulatedDielectricStatePoint),
}

impl SpectralOpticalStatePoint {
    /// Exact original source snapshot for the selected family.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Conductor(source) => source.resolved(),
            Self::TabulatedDielectric(source) => source.resolved(),
        }
    }

    /// Schema-specific identity, including all selected coordinates and claims.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        match self {
            Self::Conductor(source) => source.identity(),
            Self::TabulatedDielectric(source) => source.identity(),
        }
    }

    /// Refresh the retained family and exact selected claims together.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        match self {
            Self::Conductor(source) => source.at_query_point(card, point).map(Self::Conductor),
            Self::TabulatedDielectric(source) => source
                .at_query_point(card, point)
                .map(Self::TabulatedDielectric),
        }
    }
}

/// Resolve all 27 conductor coordinates and coefficients from one immutable
/// card at one complete state. Wavelengths are SI metres, n/k dimensionless.
pub fn resolve_spectral_conductor_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<SpectralConductorStatePoint, MaterialStatePointError> {
    let requirements = spectral_requirements(false)?;
    SpectralConductorStatePoint::try_from_material(resolve_material_state_point(
        card,
        point,
        &requirements,
        selection,
    )?)
}

/// Resolve all 27 dielectric coordinates and coefficients. Absorption requires
/// the explicit intensity-neper convention; amplitude and decibel quantities
/// cannot silently fill it. Missing data never falls back to the RGB schema.
pub fn resolve_spectral_tabulated_dielectric_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<SpectralTabulatedDielectricStatePoint, MaterialStatePointError> {
    let requirements = spectral_requirements(true)?;
    SpectralTabulatedDielectricStatePoint::try_from_material(resolve_material_state_point(
        card,
        point,
        &requirements,
        selection,
    )?)
}

fn spectral_requirements(
    dielectric: bool,
) -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let mut requirements = Vec::with_capacity(3 * SPECTRAL_OPTICAL_SAMPLE_COUNT);
    let (index_names, loss_names, loss_quantity) = if dielectric {
        (
            SPECTRAL_DIELECTRIC_IOR_PROPERTIES,
            SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
        )
    } else {
        (
            SPECTRAL_CONDUCTOR_ETA_PROPERTIES,
            SPECTRAL_CONDUCTOR_K_PROPERTIES,
            SPECTRAL_CONDUCTOR_K_QUANTITY,
        )
    };
    for (names, quantity, admissibility) in [
        (
            SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES,
            WAVELENGTH_QUANTITY,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (
            index_names,
            INDEX_QUANTITY,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (loss_names, loss_quantity, ScalarAdmissibility::NonNegative),
    ] {
        for name in names {
            requirements.push(ScalarPropertyRequirement::try_with_quantity(
                name,
                quantity,
                admissibility,
            )?);
        }
    }
    Ok(requirements)
}

fn source_wavelengths(
    resolved: &ResolvedMaterialStatePoint,
) -> Result<[f64; SPECTRAL_OPTICAL_SAMPLE_COUNT], MaterialStatePointError> {
    let mut wavelengths = [0.0; SPECTRAL_OPTICAL_SAMPLE_COUNT];
    for (index, wavelength) in wavelengths.iter_mut().enumerate() {
        *wavelength = source_value(
            resolved,
            SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[index],
            WAVELENGTH_QUANTITY,
            ScalarAdmissibility::StrictlyPositive,
        )?;
    }
    if wavelengths.windows(2).any(|pair| pair[0] >= pair[1]) {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "spectral_optical_wavelength_order",
        });
    }
    Ok(wavelengths)
}

fn source_value(
    resolved: &ResolvedMaterialStatePoint,
    name: &str,
    quantity: QuantitySpec,
    admissibility: ScalarAdmissibility,
) -> Result<f64, MaterialStatePointError> {
    let property =
        resolved
            .property(name)
            .ok_or_else(|| MaterialStatePointError::InvalidRequirement {
                property: name.to_owned(),
                reason: "finite-band optical property is missing or ambiguous",
            })?;
    if !property
        .requirement()
        .matches_key(&PropertyKey::with_quantity(name, quantity))
        || !property.is_static_at_fixed_external_conditions()
    {
        return Err(MaterialStatePointError::InvalidRequirement {
            property: name.to_owned(),
            reason: "finite-band optics require canonical scalar quantities at fixed external conditions",
        });
    }
    let value = property.value_si();
    if !admissibility.admits(value) {
        return Err(MaterialStatePointError::OutsideConsumerDomain {
            property: name.to_owned(),
            value,
            admissibility,
        });
    }
    Ok(value)
}

#[cfg(test)]
mod tests;
