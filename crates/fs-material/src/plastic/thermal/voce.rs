//! Affine sourced Voce saturation on the original thermal and native J2 owners.

use super::*;
use crate::plastic::voce::VoceJ2;
use crate::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY,
};
use fs_qty::Pressure;

#[derive(Clone, Copy, Debug)]
pub(super) struct VoceSelection {
    pub(super) saturation_temperature_slope_pa_k: f64,
    saturation_pa: f64,
    rate_per_strain: f64,
}

impl ThermoplasticSource<'_> {
    /// Select independently sourced Voce saturation on the SAME material card,
    /// caloric span and retained native plastic history. Q(T) must be globally
    /// affine and nonnegative; b must be constant and positive. Complete source
    /// lines, including internal knots and surface cuts, are admitted and pinned.
    ///
    /// The extra Helmholtz energy is Q(T) g(alpha), where
    /// g(alpha)=alpha-(1-exp(-b alpha))/b. Its entropy is -Q' g(alpha).
    /// Affine Q and constant b add no frozen-history heat-capacity derivative.
    /// A tangent-only hardening datum refuses. Selected overstress kinetics use
    /// a joint consistency solve with the same saturation energy and history.
    pub fn with_voce_hardening(
        mut self,
        selection: MaterialPropertySelection,
    ) -> Result<Self, FeedbackError> {
        if self.voce.is_some() {
            return Err(FeedbackError::Invalid("Voce hardening is already selected"));
        }
        if self.source.hardening_basis() == HardeningPropertyBasis::UniaxialTangentModulus {
            return Err(FeedbackError::Invalid(
                "Voce hardening requires an independent plastic-strain modulus, not a total-strain tangent",
            ));
        }
        let requirements = vec![
            ScalarPropertyRequirement::try_new(
                VOCE_SATURATION_INCREMENT_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::NonNegative,
            )
            .map_err(|e| FeedbackError::Source(e.into()))?,
            ScalarPropertyRequirement::try_new(
                VOCE_SATURATION_RATE_PROPERTY,
                Dims([0; 6]),
                ScalarAdmissibility::StrictlyPositive,
            )
            .map_err(|e| FeedbackError::Source(e.into()))?,
        ];
        let resolved =
            resolve_material_state_point(self.card, &self.context, &requirements, selection)
                .map_err(|e| FeedbackError::Source(e.into()))?;
        let saturation_pa = resolved
            .property(VOCE_SATURATION_INCREMENT_PROPERTY)
            .expect("the resolved Voce requirements contain Q")
            .value_si();
        let rate_per_strain = resolved
            .property(VOCE_SATURATION_RATE_PROPERTY)
            .expect("the resolved Voce requirements contain b")
            .value_si();
        let knots = self.caloric.curve().knots();
        let lower = point_at(
            &self.context,
            &self.temperature_axis,
            knots[0].temperature_k,
        )?;
        let upper = point_at(
            &self.context,
            &self.temperature_axis,
            knots[knots.len() - 1].temperature_k,
        )?;
        let mut saturation_knots = Vec::new();
        let mut additional_pins = Vec::with_capacity(requirements.len());
        for property in resolved.properties() {
            admit_rate_independent_property(property, &self.temperature_axis)?;
            let name = property.requirement().name();
            let selected = property.answer().receipt.selected;
            let claim = self
                .card
                .claims()
                .claim(selected)
                .ok_or(FeedbackError::Invalid(
                    "Voce coefficient claim left its original card",
                ))?;
            let support = self
                .card
                .claims()
                .query_line_typed(
                    &claim.key,
                    &self.temperature_axis,
                    &lower,
                    &upper,
                    ClaimSelection::Pinned(selected),
                )
                .map_err(|e| FeedbackError::Support {
                    property: name.into(),
                    reason: e.to_string(),
                })?;
            for (temperature, sample) in support.knots {
                let value = sample.evidence.value.value;
                if name == VOCE_SATURATION_INCREMENT_PROPERTY {
                    if !value.is_finite() || value < 0.0 {
                        return Err(FeedbackError::Invalid(
                            "Voce saturation must remain nonnegative over the complete thermal span",
                        ));
                    }
                    saturation_knots.push((temperature, value));
                } else if !positive(value) || value != rate_per_strain {
                    return Err(FeedbackError::Invalid(
                        "Voce b must remain positive and constant over the complete thermal span",
                    ));
                }
            }
            additional_pins.push((name.to_owned(), selected));
        }
        let saturation_temperature_slope_pa_k = affine_temperature_slope(
            &saturation_knots,
            self.reference_temperature_k,
            saturation_pa,
        )
        .map_err(FeedbackError::Invalid)?;
        let MaterialPropertySelection::PinnedByProperty(pins) = &mut self.selection else {
            return Err(FeedbackError::Invalid(
                "thermoplastic coefficients lost original claim pins",
            ));
        };
        pins.extend(additional_pins);
        self.requirements.extend(requirements);
        self.voce = Some(VoceSelection {
            saturation_temperature_slope_pa_k,
            saturation_pa,
            rate_per_strain,
        });
        self.source = self.resolve_temperature(self.reference_temperature_k)?;
        // The original native constructor admits the combined return modulus
        // at both affine extrema, including an exactly zero-Q endpoint.
        self.resolve_temperature(knots[0].temperature_k)?;
        self.resolve_temperature(knots[knots.len() - 1].temperature_k)?;
        Ok(self)
    }

    /// Whether nonlinear energetic saturation was explicitly selected.
    pub const fn has_voce_hardening(&self) -> bool {
        self.voce.is_some()
    }

    pub(super) fn voce_law<E>(
        &self,
        source: &ResolvedJ2Plasticity,
    ) -> Result<Option<VoceJ2>, FeedbackError<E>> {
        let Some(selection) = self.voce else {
            return Ok(None);
        };
        let value = |name| {
            source
                .material()
                .property(name)
                .map(|p| p.value_si())
                .ok_or(FeedbackError::Invalid(
                    "resolved source omitted a Voce coefficient",
                ))
        };
        let saturation = value(VOCE_SATURATION_INCREMENT_PROPERTY)?;
        let rate = value(VOCE_SATURATION_RATE_PROPERTY)?;
        if rate != selection.rate_per_strain
            || selection.saturation_temperature_slope_pa_k == 0.0
                && saturation != selection.saturation_pa
        {
            return Err(FeedbackError::Invalid(
                "pinned coefficients left the admitted affine Voce free-energy model",
            ));
        }
        VoceJ2::new(source.law().clone(), saturation, rate)
            .map(Some)
            .map_err(|e| FeedbackError::Source(e.into()))
    }
}

impl CaloricPoint<'_> {
    /// Select sourced Voce hardening before loading or advancing the specimen.
    /// Its physical heat and reversible exchange use the original caloric owner.
    pub fn with_voce_hardening(
        mut self,
        selection: MaterialPropertySelection,
    ) -> Result<Self, FeedbackError> {
        if self.time_s != 0.0 || self.strain != [0.0; 6] || self.history.alpha != 0.0 {
            return Err(FeedbackError::Invalid(
                "select Voce hardening before advancing the specimen",
            ));
        }
        self.binding = self.binding.with_voce_hardening(selection)?;
        self.source = self
            .binding
            .resolve_temperature(self.state.temperature_k())?;
        Ok(self)
    }
}
