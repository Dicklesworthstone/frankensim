//! Source-backed Sellmeier coefficients for the existing visible dielectric path.
//!
//! Coefficients belong to one immutable material card and one physical state.
//! Wavelengths in the source equation are vacuum metres:
//! `n² = 1 + sum B_i lambda² / (lambda² - C_i)`.
//! The selected source must express absolute index relative to vacuum; catalog
//! indices relative to ambient air require a source-side conversion first.
//! This resolver checks source support, units and scalar domains. The optical
//! consumer additionally admits its transparent wavelength band and rejects
//! resonances or nonphysical indices there. No optical constants are inferred
//! from material names.

use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::Dims;

use super::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, VisibleDielectricAbsorption,
    append_visible_dielectric_absorption_requirements, resolve_material_state_point,
    visible_dielectric_absorption, visible_dielectric_uses_spectral_absorption,
};

/// Dimensionless oscillator strengths for `n² = 1 + sum B_i lambda²/(lambda²-C_i)`.
pub const VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES: [&str; 3] = [
    "optical_sellmeier_b1",
    "optical_sellmeier_b2",
    "optical_sellmeier_b3",
];
/// Squared vacuum resonance wavelengths in SI m², not micrometres squared.
pub const VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES: [&str; 3] = [
    "optical_sellmeier_c1_m2",
    "optical_sellmeier_c2_m2",
    "optical_sellmeier_c3_m2",
];

/// A homogeneous transparent dielectric source with three Sellmeier terms.
///
/// Nonnegative oscillator strengths and squared resonance wavelengths describe
/// the lossless normal-dispersion model between resonances. Absorption remains
/// the explicitly declared spectral or RGB Beer-Lambert response; these
/// coefficients alone do not model resonant loss or anomalous dispersion.
#[derive(Clone, Debug, PartialEq)]
pub struct VisibleSellmeierStatePoint {
    resolved: ResolvedMaterialStatePoint,
    strengths: [f64; 3],
    resonances_m2: [f64; 3],
    absorption: VisibleDielectricAbsorption,
}

impl VisibleSellmeierStatePoint {
    /// Complete source card, query coordinates and property-use receipts.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Dimensionless strengths and squared resonance wavelengths [m²].
    #[must_use]
    pub const fn sellmeier_coefficients_si(&self) -> ([f64; 3], [f64; 3]) {
        (self.strengths, self.resonances_m2)
    }

    /// Absorption in its original source parameterization.
    #[must_use]
    pub const fn absorption(&self) -> VisibleDielectricAbsorption {
        self.absorption
    }

    /// Refresh every oscillator and absorption coefficient from its original
    /// selected claim at a new complete physical condition. The original
    /// spectral or RGB absorption schema and scalar bounds remain binding.
    /// The renderer must re-admit the refreshed transparent wavelength band;
    /// source support alone does not exclude a resonance moving into that band.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        visible_sellmeier_from_resolved(
            self.resolved.at_query_point(card, point)?,
            matches!(
                self.absorption,
                VisibleDielectricAbsorption::Spectral { .. }
            ),
        )
    }
}

/// Resolve all six Sellmeier coefficients and the complete absorption schema.
///
/// Every property must be supported by the same card at the exact supplied
/// state, with the original claim-selection policy. Missing coefficients,
/// wrong dimensions, extrapolation and gain refuse atomically. Zero-strength
/// terms must still be declared explicitly; no coefficient is defaulted.
/// Absorption uses the same typed spectral or RGB source-selection rules as
/// [`super::resolve_visible_dielectric_state_point`], including ambiguity and
/// refusal of incomplete spectral data without an RGB fallback.
pub fn resolve_visible_sellmeier_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<VisibleSellmeierStatePoint, MaterialStatePointError> {
    let spectral = visible_dielectric_uses_spectral_absorption(card, &selection)?;
    let mut requirements = Vec::with_capacity(if spectral { 15 } else { 10 });
    for name in VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES {
        requirements.push(ScalarPropertyRequirement::try_new(
            name,
            Dims::NONE,
            ScalarAdmissibility::NonNegative,
        )?);
    }
    for name in VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES {
        requirements.push(ScalarPropertyRequirement::try_new(
            name,
            Dims([2, 0, 0, 0, 0, 0]),
            ScalarAdmissibility::NonNegative,
        )?);
    }
    append_visible_dielectric_absorption_requirements(&mut requirements, spectral)?;
    let resolved = resolve_material_state_point(card, point, &requirements, selection)?;
    visible_sellmeier_from_resolved(resolved, spectral)
}

fn visible_sellmeier_from_resolved(
    resolved: ResolvedMaterialStatePoint,
    spectral: bool,
) -> Result<VisibleSellmeierStatePoint, MaterialStatePointError> {
    let value = |name: &str| {
        resolved
            .property(name)
            .expect("canonical Sellmeier requirement was resolved")
            .value_si()
    };
    let strengths = VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES.map(value);
    let resonances_m2 = VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES.map(value);
    let absorption = visible_dielectric_absorption(&resolved, spectral)?;
    Ok(VisibleSellmeierStatePoint {
        resolved,
        strengths,
        resonances_m2,
        absorption,
    })
}
