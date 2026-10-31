//! Source-derived two-temperature Hertz compliance and its thermal derivatives.
//!
//! Each original solid keeps its own temperature, complete fixed external
//! context, immutable card and selected claims. E_a(T_a) and E_b(T_b) are
//! affine on admitted source intervals; each Poisson ratio and yield limit is
//! constant there. The ORIGINAL reduced-modulus owner evaluates
//! E* = [(1-nu_a^2)/E_a + (1-nu_b^2)/E_b]^-1. E* is generally NOT affine.
//!
//! Derivatives are at fixed indentation/geometry. If the native normal law has
//! psi=G(delta)*E*, its surface entropies are -G*dE*/dTi and its internal
//! contribution is psi + sum(Ti*si). The Hessian supplies the reciprocal heat
//! capacity. This module supplies no contact integrator, geometry update,
//! measured calibration, plasticity, phase change or expansion law.

use crate::state_point::{
    DENSITY_PROPERTY, IsotropicSolidStatePoint, MaterialPropertySelection, MaterialStatePointError,
    POISSON_RATIO_PROPERTY, ResolvedReducedModulus, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    affine_temperature_slope, resolve_isotropic_solid_state_point,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

#[derive(Debug)]
pub enum ThermalHertzError {
    Material(MaterialStatePointError),
    Invalid(&'static str),
}
impl From<MaterialStatePointError> for ThermalHertzError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for ThermalHertzError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-derived thermal Hertz compliance: {self:?}")
    }
}
impl std::error::Error for ThermalHertzError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material(error) => Some(error),
            Self::Invalid(_) => None,
        }
    }
}
fn finite(value: f64) -> Result<f64, ThermalHertzError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ThermalHertzError::Invalid(
            "unrepresentable thermal contact derivative",
        ))
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Surface {
    card: MaterialCard,
    reference: IsotropicSolidStatePoint,
    temperature_range_k: [f64; 2],
    reference_temperature_k: f64,
    youngs_slope_pa_k: f64,
}
impl Surface {
    fn new(
        card: &MaterialCard,
        reference: &IsotropicSolidStatePoint,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, ThermalHertzError> {
        let invalid = ThermalHertzError::Invalid;
        if reference.elastic_state().derivation().is_some() {
            return Err(invalid("affine thermal contact requires its original canonical Young/Poisson source law"));
        }
        let [low, high] = temperature_range_k;
        if !low.is_finite() || !high.is_finite() || low <= 0.0 || low >= high {
            return Err(invalid(
                "positive increasing source temperature interval required",
            ));
        }
        let original = reference.resolved();
        // This original operation pins the card, source claims and query schema.
        original.at_query_point(card, &original.query())?;
        let reference_temperature_k = *original
            .query()
            .axes()
            .get("T")
            .ok_or(invalid("canonical T coordinate is required"))?;
        if !(low..=high).contains(&reference_temperature_k) {
            return Err(invalid(
                "initial temperature leaves the admitted source interval",
            ));
        }
        let temperature_quantity = *original.axis_quantities().get("T").ok_or(invalid(
            "temperature needs explicit kelvin dimensions or absolute semantics",
        ))?;
        if temperature_quantity != QuantitySpec::dimensional(Temperature::DIMS)
            && temperature_quantity
                != QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                ))
        {
            return Err(invalid("contact temperature must be absolute kelvin"));
        }
        let mut youngs_slope_pa_k = 0.0;
        for name in [
            DENSITY_PROPERTY,
            YOUNG_MODULUS_PROPERTY,
            POISSON_RATIO_PROPERTY,
            YIELD_STRESS_PROPERTY,
        ] {
            let property = original.property(name).ok_or(invalid(
                "original canonical solid property is missing or ambiguous",
            ))?;
            if !property.has_only_owned_and_fixed_external_conditions("T", temperature_quantity) {
                return Err(invalid(
                    "thermal contact requires one owned temperature and fixed external conditions",
                ));
            }
            let claim = card
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or(invalid("original solid source claim is absent"))?;
            if claim
                .value
                .axes()
                .any(|(axis, dims)| axis == "T" && dims != Temperature::DIMS)
            {
                return Err(invalid(
                    "source temperature axis has incompatible dimensions",
                ));
            }
            let point = |t| {
                original
                    .query()
                    .with("T", t)
                    .map_err(|_| invalid("thermal source query is unrepresentable"))
            };
            let line = card
                .claims()
                .query_line_typed(
                    &claim.key,
                    "T",
                    &point(low)?,
                    &point(high)?,
                    ClaimSelection::Pinned(property.answer().receipt.selected),
                )
                .map_err(|source| MaterialStatePointError::Support {
                    property: name.to_owned(),
                    source: Box::new(source),
                })?;
            let knots: Vec<_> = line
                .knots
                .into_iter()
                .map(|(t, answer)| (t, answer.evidence.value.value))
                .collect();
            if knots.first().map(|k| k.0) != Some(low)
                || knots.last().map(|k| k.0) != Some(high)
                || knots.iter().any(|(_, v)| {
                    !v.is_finite()
                        || if name == POISSON_RATIO_PROPERTY {
                            *v <= -1.0 || *v >= 0.5
                        } else {
                            *v <= 0.0
                        }
                })
            {
                return Err(invalid(
                    "solid source must continuously cover the complete positive thermal interval",
                ));
            }
            if name == YOUNG_MODULUS_PROPERTY {
                youngs_slope_pa_k =
                    affine_temperature_slope(&knots, reference_temperature_k, property.value_si())
                        .map_err(invalid)?;
            } else if name == POISSON_RATIO_PROPERTY || name == YIELD_STRESS_PROPERTY {
                // Constant nu keeps the native Hertz compliance differentiable
                // with only the admitted E slopes. The unchanged yield limit
                // remains valid for the original elastic applicability gate.
                if knots
                    .iter()
                    .any(|(_, v)| v.to_bits() != property.value_si().to_bits())
                {
                    return Err(invalid(
                        "this thermal Hertz rung needs constant Poisson ratio and yield limit",
                    ));
                }
            }
        }
        let value = Self {
            card: card.clone(),
            reference: reference.clone(),
            temperature_range_k,
            reference_temperature_k,
            youngs_slope_pa_k,
        };
        for t in [low, reference_temperature_k, high] {
            value.at_temperature(t)?;
        }
        Ok(value)
    }

    fn at_temperature(&self, t: f64) -> Result<IsotropicSolidStatePoint, ThermalHertzError> {
        if !t.is_finite()
            || !(self.temperature_range_k[0]..=self.temperature_range_k[1]).contains(&t)
        {
            return Err(ThermalHertzError::Invalid(
                "temperature leaves the pinned solid source interval",
            ));
        }
        let point = self
            .reference
            .resolved()
            .query()
            .with("T", t)
            .map_err(|_| ThermalHertzError::Invalid("thermal source query is unrepresentable"))?;
        let selection = MaterialPropertySelection::PinnedByProperty(
            self.reference
                .resolved()
                .properties()
                .iter()
                .map(|p| {
                    (
                        p.requirement().name().to_owned(),
                        p.answer().receipt.selected,
                    )
                })
                .collect(),
        );
        let result = resolve_isotropic_solid_state_point(&self.card, &point, selection)?;
        let expected = finite(self.youngs_slope_pa_k.mul_add(
            t - self.reference_temperature_k,
            self.reference.young_modulus_pa(),
        ))?;
        if (expected - result.young_modulus_pa()).abs()
            > 512.0 * f64::EPSILON * expected.abs().max(result.young_modulus_pa().abs())
            || result.poisson_ratio().to_bits() != self.reference.poisson_ratio().to_bits()
            || result.yield_stress_pa().to_bits() != self.reference.yield_stress_pa().to_bits()
        {
            return Err(ThermalHertzError::Invalid(
                "refreshed solid differs from its admitted affine/fixed sources",
            ));
        }
        Ok(result)
    }
}

/// Original two-solid compliance on a complete rectangular temperature domain.
/// The source laws are immutable; only actual surface temperatures evolve.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalHertzPair {
    surfaces: std::sync::Arc<[Surface; 2]>,
    identity: ContentHash,
}

/// Native queried states and exact analytical derivatives of the ORIGINAL
/// reduced modulus, in ordered surface-A/surface-B temperature coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermalHertz {
    pub surfaces: [IsotropicSolidStatePoint; 2],
    pub reduced_modulus: ResolvedReducedModulus,
    pub gradient_pa_k: [f64; 2],
    pub hessian_pa_k2: [[f64; 2]; 2],
}
impl ThermalHertzPair {
    pub fn new(
        cards: [&MaterialCard; 2],
        references: [&IsotropicSolidStatePoint; 2],
        temperature_ranges_k: [[f64; 2]; 2],
    ) -> Result<Self, ThermalHertzError> {
        let surfaces = [
            Surface::new(cards[0], references[0], temperature_ranges_k[0])?,
            Surface::new(cards[1], references[1], temperature_ranges_k[1])?,
        ];
        let mut hash = DomainHasher::new("org.frankensim.fs-material.thermal-hertz-pair.v1");
        for surface in &surfaces {
            hash.update(surface.card.content_hash().as_bytes());
            hash.update(surface.reference.resolved().identity().as_bytes());
            for value in surface
                .temperature_range_k
                .into_iter()
                .chain([surface.youngs_slope_pa_k])
            {
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
        let result = Self {
            surfaces: std::sync::Arc::new(surfaces),
            identity: hash.finalize(),
        };
        for a in temperature_ranges_k[0] {
            for b in temperature_ranges_k[1] {
                result.at_temperatures([a, b])?;
            }
        }
        Ok(result)
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub fn card_identities(&self) -> [ContentHash; 2] {
        self.surfaces
            .each_ref()
            .map(|surface| surface.card.content_hash())
    }
    pub fn reference_temperatures_k(&self) -> [f64; 2] {
        self.surfaces
            .each_ref()
            .map(|surface| surface.reference_temperature_k)
    }
    pub fn temperature_ranges_k(&self) -> [[f64; 2]; 2] {
        self.surfaces
            .each_ref()
            .map(|surface| surface.temperature_range_k)
    }
    pub fn limiting_yield_stress_pa(&self) -> f64 {
        self.surfaces[0]
            .reference
            .yield_stress_pa()
            .min(self.surfaces[1].reference.yield_stress_pa())
    }

    pub fn at_temperatures(
        &self,
        temperatures_k: [f64; 2],
    ) -> Result<ResolvedThermalHertz, ThermalHertzError> {
        let surfaces = [
            self.surfaces[0].at_temperature(temperatures_k[0])?,
            self.surfaces[1].at_temperature(temperatures_k[1])?,
        ];
        let reduced_modulus = surfaces[0].reduced_modulus_against(&surfaces[1])?;
        let e = reduced_modulus.value_pa;
        let mut weights = [0.0; 2];
        let mut ratios = [0.0; 2];
        let mut gradient_pa_k = [0.0; 2];
        for i in 0..2 {
            let bulk = surfaces[i].young_modulus_pa();
            weights[i] = finite(
                ((1.0 - surfaces[i].poisson_ratio() * surfaces[i].poisson_ratio()) / bulk) * e,
            )?;
            ratios[i] = finite(self.surfaces[i].youngs_slope_pa_k / bulk)?;
            gradient_pa_k[i] = finite(finite(e * weights[i])? * ratios[i])?;
            if self.surfaces[i].youngs_slope_pa_k != 0.0
                && (ratios[i] == 0.0 || gradient_pa_k[i] == 0.0)
            {
                return Err(ThermalHertzError::Invalid(
                    "nonzero source derivative underflows",
                ));
            }
        }
        // Concavity of harmonic compliance gives a negative-semidefinite
        // Hessian. The factored form avoids subtracting nearly equal terms.
        let scale = finite(finite(2.0 * e * weights[0])? * weights[1])?;
        let cross = finite(finite(scale * ratios[0])? * ratios[1])?;
        let hessian_pa_k2 = [
            [finite(-finite(scale * ratios[0])? * ratios[0])?, cross],
            [cross, finite(-finite(scale * ratios[1])? * ratios[1])?],
        ];
        Ok(ResolvedThermalHertz {
            surfaces,
            reduced_modulus,
            gradient_pa_k,
            hessian_pa_k2,
        })
    }
}

#[cfg(test)]
mod tests;
