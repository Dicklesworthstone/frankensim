//! Original-source temperature derivatives of a native Neo-Hookean model.
//!
//! The caller qualifies E(T) and nu(T) as the infinitesimal moduli of this
//! thermally relaxed material. The same-card, same-context sources must match
//! the selected native mu and lambda at the actual temperature. Their admitted
//! affine paths supply exact rational coefficient derivatives; this does not
//! refit a model, reinterpret a prestressed tangent modulus, or change its
//! current energy/stress/tangent. Source receipts remain attached.

use super::{ElasticTemperatureDerivatives, FiniteThermoelasticState, ThermalReferenceJet};
use crate::MaterialError;
use crate::hyper::cards::ResolvedHyperelastic;
use crate::hyper::{Hyperelastic, HyperelasticModel, neo_hookean_energy};
use crate::state_point::{
    MaterialStatePointError, POISSON_RATIO_PROPERTY, ResolvedMaterialStatePoint,
    YOUNG_MODULUS_PROPERTY, affine_temperature_slope,
};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyLineSlice};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

/// Missing, incompatible or nondifferentiable original coefficient evidence.
#[derive(Clone, Debug, PartialEq)]
pub enum NeoHookeanThermalSourceError {
    Invalid(&'static str),
    Source(MaterialStatePointError),
    Material(MaterialError),
}
impl core::fmt::Display for NeoHookeanThermalSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Neo-Hookean thermal source: {self:?}")
    }
}
impl std::error::Error for NeoHookeanThermalSourceError {}
impl From<MaterialStatePointError> for NeoHookeanThermalSourceError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Source(error)
    }
}
impl From<MaterialError> for NeoHookeanThermalSourceError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
type Result<T> = core::result::Result<T, NeoHookeanThermalSourceError>;

fn invalid(message: &'static str) -> NeoHookeanThermalSourceError {
    NeoHookeanThermalSourceError::Invalid(message)
}
fn same_coefficient(left: f64, right: f64) -> bool {
    left.is_finite()
        && right.is_finite()
        && (left == right
            || (left - right).abs() <= 128.0 * f64::EPSILON * left.abs().max(right.abs()))
}

struct CoefficientJet {
    values: [f64; 2],
    first: [f64; 2],
    second: [f64; 2],
}

fn coefficient_jet(values: [f64; 2], slopes: [f64; 2]) -> Result<CoefficientJet> {
    let [youngs, poisson] = values;
    let [youngs_t, poisson_t] = slopes;
    if !youngs.is_finite() || youngs <= 0.0 || !(-1.0..0.5).contains(&poisson) {
        return Err(invalid(
            "original E/nu coefficients are outside their elastic domain",
        ));
    }
    let shear = youngs / (2.0 * (1.0 + poisson));
    let ratio = poisson / (1.0 - 2.0 * poisson);
    let lame = 2.0 * shear * ratio;
    let (first, second) = if slopes == [0.0; 2] {
        ([0.0; 2], [0.0; 2])
    } else {
        let poisson_shear_rate = poisson_t / (1.0 + poisson);
        let shear_t = youngs_t / (2.0 * (1.0 + poisson)) - shear * poisson_shear_rate;
        let shear_tt = -2.0 * poisson_shear_rate * shear_t;
        let poisson_bulk_rate = poisson_t / (1.0 - 2.0 * poisson);
        let ratio_t = poisson_bulk_rate / (1.0 - 2.0 * poisson);
        let ratio_tt = 4.0 * poisson_bulk_rate * ratio_t;
        let lame_t = 2.0 * (shear_t * ratio + shear * ratio_t);
        let lame_tt = 2.0 * (shear_tt * ratio + 2.0 * shear_t * ratio_t + shear * ratio_tt);
        ([shear_t, lame_t], [shear_tt, lame_tt])
    };
    if [shear, lame]
        .iter()
        .chain(&first)
        .chain(&second)
        .any(|value| !value.is_finite())
        || (slopes != [0.0; 2] && first == [0.0; 2])
    {
        return Err(invalid(
            "original coefficient temperature derivatives are unrepresentable",
        ));
    }
    Ok(CoefficientJet {
        values: [shear, lame],
        first,
        second,
    })
}

/// One source/model snapshot with exact coefficient derivatives at its actual T.
///
/// E and nu are independently affine over the explicit source interval. Their
/// derived Lamé coefficients generally have nonzero second derivatives. This
/// The original model remains the anchor snapshot. Explicit temperature
/// evaluation derives a new native response from these same pinned source
/// paths; it does not mutate or relabel the anchor model card.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedNeoHookeanThermalCoefficients {
    model: ResolvedHyperelastic,
    material: ResolvedMaterialStatePoint,
    source_lines: [PropertyLineSlice; 2],
    temperature_range_k: [f64; 2],
    temperature_k: f64,
    source_values: [f64; 2],
    source_slopes: [f64; 2],
    first_derivatives_pa_k: [f64; 2],
    second_derivatives_pa_k2: [f64; 2],
}

impl ResolvedNeoHookeanThermalCoefficients {
    /// Retain the original model, point and complete pinned E/nu source lines.
    ///
    /// Both source paths must be smooth and cover the whole positive interval,
    /// with the actual temperature strictly inside it. Every nonthermal
    /// coordinate stays fixed. The native Neo-Hookean model must belong to the
    /// same card and agree with the derived current coefficients to roundoff;
    /// a discrepant measurement is refused, never averaged or silently fitted.
    pub fn try_new(
        card: &MaterialCard,
        model: &ResolvedHyperelastic,
        material: &ResolvedMaterialStatePoint,
        temperature_range_k: [f64; 2],
    ) -> Result<Self> {
        let HyperelasticModel::NeoHookean { mu, lambda } = model.law().model else {
            return Err(invalid(
                "thermal E/nu derivatives require the original Neo-Hookean family",
            ));
        };
        let point = material.query();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let temperature_k = *point
            .axes()
            .get("T")
            .ok_or_else(|| invalid("the original absolute temperature coordinate is missing"))?;
        let [low, high] = temperature_range_k;
        if material.card_identity() != card.content_hash()
            || model.context() != &point
            || !card
                .models()
                .iter()
                .any(|source| source == model.model_card())
            || point.axis_quantities().get("T") != Some(&absolute)
            || !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || !(low < temperature_k && temperature_k < high)
        {
            return Err(invalid(
                "same-card model and typed source context require a two-sided positive temperature interval",
            ));
        }
        material.at_query_point(card, &point)?;
        let lower = point
            .clone()
            .with("T", low)
            .map_err(|_| invalid("the original lower temperature point is invalid"))?;
        let upper = point
            .clone()
            .with("T", high)
            .map_err(|_| invalid("the original upper temperature point is invalid"))?;
        if [&lower, &upper].iter().any(|point| {
            !model
                .model_card()
                .validity
                .contains_typed(point.axes(), point.axis_quantities())
        }) {
            return Err(invalid(
                "the derivative interval leaves the original model domain",
            ));
        }
        let mut values = [0.0; 2];
        let mut slopes = [0.0; 2];
        let mut lines = Vec::with_capacity(2);
        for (index, (name, dimensions)) in [
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS),
            (POISSON_RATIO_PROPERTY, Dims::NONE),
        ]
        .into_iter()
        .enumerate()
        {
            let property = material
                .property(name)
                .ok_or_else(|| invalid("an original E or nu source is missing or ambiguous"))?;
            if property.requirement().quantity() != QuantitySpec::dimensional(dimensions)
                || !property.has_only_fixed_external_conditions(Some("T"), false)
            {
                return Err(invalid(
                    "original static E/nu meanings and fixed external coordinates are required",
                ));
            }
            let selected = property.answer().receipt.selected;
            let claim = card
                .claims()
                .claim(selected)
                .ok_or_else(|| invalid("an original coefficient claim is missing"))?;
            if claim
                .value
                .axes()
                .any(|(axis, dims)| axis == "T" && dims != Temperature::DIMS)
            {
                return Err(invalid(
                    "the source temperature axis has incompatible dimensions",
                ));
            }
            let line = card
                .claims()
                .query_line_typed(
                    &claim.key,
                    "T",
                    &lower,
                    &upper,
                    ClaimSelection::Pinned(selected),
                )
                .map_err(|source| MaterialStatePointError::Support {
                    property: name.to_owned(),
                    source: Box::new(source),
                })?;
            let knots: Vec<_> = line
                .knots
                .iter()
                .map(|(t, answer)| (*t, answer.evidence.value.value))
                .collect();
            if knots.first().map(|knot| knot.0) != Some(low)
                || knots.last().map(|knot| knot.0) != Some(high)
                || knots.iter().any(|(_, value)| {
                    !value.is_finite()
                        || if index == 0 {
                            *value <= 0.0
                        } else {
                            *value <= -1.0 || *value >= 0.5
                        }
                })
            {
                return Err(invalid(
                    "the complete affine E/nu source interval must remain physically admissible",
                ));
            }
            values[index] = property.value_si();
            // Use the same exact source-line and slope admission as the native
            // small-strain thermodynamic owners, without finite differences.
            let (reference_t, reference_value) = if index == 0 {
                (temperature_k, values[index])
            } else {
                knots[0]
            };
            slopes[index] =
                affine_temperature_slope(&knots, reference_t, reference_value).map_err(invalid)?;
            let expected = slopes[index].mul_add(temperature_k - knots[0].0, knots[0].1);
            let scale = knots.iter().fold(
                expected.abs().max(values[index].abs()),
                |scale, (_, value)| scale.max(value.abs()),
            );
            if !expected.is_finite()
                || (expected - values[index]).abs() > 128.0 * f64::EPSILON * scale
                || (slopes[index] == 0.0 && knots[0].1 != knots[knots.len() - 1].1)
            {
                return Err(invalid(
                    "the current coefficient or derivative differs from its original source path",
                ));
            }
            lines.push(line);
        }
        let jet = coefficient_jet(values, slopes)?;
        if !same_coefficient(jet.values[0], mu) || !same_coefficient(jet.values[1], lambda) {
            return Err(invalid(
                "original E/nu sources disagree with the selected native mu/lambda",
            ));
        }
        let mut lines = lines.into_iter();
        Ok(Self {
            model: model.clone(),
            material: material.clone(),
            source_lines: [
                lines.next().expect("E source"),
                lines.next().expect("nu source"),
            ],
            temperature_range_k,
            temperature_k,
            source_values: values,
            source_slopes: slopes,
            first_derivatives_pa_k: jet.first,
            second_derivatives_pa_k2: jet.second,
        })
    }

    pub const fn model(&self) -> &ResolvedHyperelastic {
        &self.model
    }
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Complete original E and nu source-line answers and receipts, in that order.
    pub const fn source_lines(&self) -> &[PropertyLineSlice; 2] {
        &self.source_lines
    }
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }
    /// Exact partial derivatives [dmu/dT, dlambda/dT] of the original sources.
    pub const fn first_derivatives_pa_k(&self) -> [f64; 2] {
        self.first_derivatives_pa_k
    }
    /// Exact partial derivatives [d²mu/dT², d²lambda/dT²].
    pub const fn second_derivatives_pa_k2(&self) -> [f64; 2] {
        self.second_derivatives_pa_k2
    }

    /// Original native finite potential with its complete thermal chain rule.
    /// The supplied natural-map jet retains its independent source convention.
    /// Caloric storage is added separately by the existing adiabatic consumer.
    pub fn thermoelastic_state(
        &self,
        total_deformation_gradient: &[f64; 9],
        reference_jet: ThermalReferenceJet,
    ) -> Result<FiniteThermoelasticState> {
        thermal_state(
            self.model.law(),
            self.temperature_k,
            self.first_derivatives_pa_k,
            self.second_derivatives_pa_k2,
            total_deformation_gradient,
            reference_jet,
        )
    }

    /// Evaluate a new native response on the original affine E(T)/nu(T) paths.
    /// The caller supplies the natural-map jet at this SAME temperature. All
    /// original source receipts, fixed conditions and the anchor model remain
    /// retained on this owner; no property re-selection or numerical derivative
    /// is performed. The original elastic stretch envelope still gates Fe.
    /// Interval endpoints use the admitted path's one-sided thermal derivative.
    pub fn thermoelastic_state_at_temperature(
        &self,
        temperature_k: f64,
        total_deformation_gradient: &[f64; 9],
        reference_jet: ThermalReferenceJet,
    ) -> Result<FiniteThermoelasticState> {
        if !temperature_k.is_finite()
            || temperature_k < self.temperature_range_k[0]
            || temperature_k > self.temperature_range_k[1]
        {
            return Err(invalid(
                "temperature leaves the original coefficient source interval",
            ));
        }
        if temperature_k == self.temperature_k {
            return self.thermoelastic_state(total_deformation_gradient, reference_jet);
        }
        if self.source_slopes == [0.0; 2] {
            return thermal_state(
                self.model.law(),
                temperature_k,
                [0.0; 2],
                [0.0; 2],
                total_deformation_gradient,
                reference_jet,
            );
        }
        let values = core::array::from_fn(|i| {
            self.source_slopes[i].mul_add(temperature_k - self.temperature_k, self.source_values[i])
        });
        let jet = coefficient_jet(values, self.source_slopes)?;
        let law = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: jet.values[0],
                lambda: jet.values[1],
            },
            self.model.law().stretch_limit,
        )?;
        thermal_state(
            &law,
            temperature_k,
            jet.first,
            jet.second,
            total_deformation_gradient,
            reference_jet,
        )
    }
}

fn thermal_state(
    law: &Hyperelastic,
    temperature_k: f64,
    first: [f64; 2],
    second: [f64; 2],
    total_deformation_gradient: &[f64; 9],
    reference_jet: ThermalReferenceJet,
) -> Result<FiniteThermoelasticState> {
    if first == [0.0; 2] && second == [0.0; 2] {
        return Ok(FiniteThermoelasticState::try_new(
            law,
            total_deformation_gradient,
            temperature_k,
            reference_jet,
        )?);
    }
    Ok(
        FiniteThermoelasticState::with_elastic_temperature_derivatives(
            law,
            total_deformation_gradient,
            temperature_k,
            reference_jet,
            |elastic_f| {
                let [mu_t, lambda_t] = first;
                let [mu_tt, lambda_tt] = second;
                // The original native potential is linear in mu and lambda.
                // Signed derivative coefficients are not substitute physical
                // laws; the original positive-modulus owner admits actual Fe.
                let (energy_j_m3_k, piola_pa_k) =
                    fs_ad::gradient::<9>(*elastic_f, |f| neo_hookean_energy(&f, mu_t, lambda_t));
                Ok(ElasticTemperatureDerivatives {
                    energy_j_m3_k,
                    energy_j_m3_k2: neo_hookean_energy(elastic_f, mu_tt, lambda_tt),
                    piola_pa_k,
                })
            },
        )?,
    )
}
