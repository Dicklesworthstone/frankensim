//! Original-source selection for linear time recovery of each cyclic branch.
//! Rates are independent kinetic data, not inferred from dynamic gamma,
//! elasticity, melting point, viscosity or a temperature activation heuristic.

use super::*;
use crate::plastic::chaboche::recovery::STATIC_RECOVERY_RATE_DIMS;

impl<'a, const N: usize> ThermalChabocheSource<'a, N> {
    /// Select one original nonnegative rate [s^-1] for EVERY ordered branch.
    /// Each coefficient is the complete r in dX/dt=-r X. Source curves may
    /// vary with temperature; their full support must cover the caloric span.
    /// Fixed external context, source pins and branch identity are preserved.
    /// Actual dt is mandatory in every subsequent constitutive proposal.
    pub fn with_static_recovery(
        mut self,
        selection: MaterialPropertySelection,
        properties: [String; N],
    ) -> Result<Self, FeedbackError> {
        if self.static_recovery_properties.is_some() {
            return Err(FeedbackError::Invalid(
                "selected static recovery cannot be replaced underneath a cyclic family",
            ));
        }
        let mut requirements = Vec::with_capacity(N);
        for (i, name) in properties.iter().enumerate() {
            if self
                .base
                .requirements
                .iter()
                .any(|entry| entry.name() == name)
                || properties[..i].contains(name)
            {
                return Err(FeedbackError::Invalid(
                    "static-recovery names must be distinct and cannot alias existing source roles",
                ));
            }
            requirements.push(
                ScalarPropertyRequirement::try_new(
                    name.clone(),
                    STATIC_RECOVERY_RATE_DIMS,
                    ScalarAdmissibility::NonNegative,
                )
                .map_err(|error| FeedbackError::Source(error.into()))?,
            );
        }
        let resolved = resolve_material_state_point(
            self.base.card,
            &self.base.context,
            &requirements,
            selection,
        )
        .map_err(|error| FeedbackError::Source(error.into()))?;
        let knots = self.base.caloric.curve().knots();
        let lower = point_at(
            &self.base.context,
            &self.base.temperature_axis,
            knots[0].temperature_k,
        )?;
        let upper = point_at(
            &self.base.context,
            &self.base.temperature_axis,
            knots[knots.len() - 1].temperature_k,
        )?;
        let mut additional_pins = Vec::with_capacity(N);
        for name in &properties {
            let property = resolved.property(name).ok_or(FeedbackError::Invalid(
                "resolved source omitted a static-recovery rate",
            ))?;
            admit_rate_independent_property(property, &self.base.temperature_axis)?;
            let selected = property.answer().receipt.selected;
            let claim = self
                .base
                .card
                .claims()
                .claim(selected)
                .ok_or(FeedbackError::Invalid(
                    "static recovery left its original material card",
                ))?;
            let support = self
                .base
                .card
                .claims()
                .query_line_typed(
                    &claim.key,
                    &self.base.temperature_axis,
                    &lower,
                    &upper,
                    ClaimSelection::Pinned(selected),
                )
                .map_err(|error| FeedbackError::Support {
                    property: name.clone(),
                    reason: error.to_string(),
                })?;
            if support.knots.iter().any(|(_, sample)| {
                let rate = sample.evidence.value.value;
                !rate.is_finite() || rate < 0.0
            }) {
                return Err(FeedbackError::Invalid(
                    "static-recovery rates must remain finite and nonnegative throughout source support",
                ));
            }
            additional_pins.push((name.clone(), selected));
        }
        let MaterialPropertySelection::PinnedByProperty(pins) = &mut self.base.selection else {
            return Err(FeedbackError::Invalid(
                "thermal cyclic source lost its original claim pins",
            ));
        };
        pins.extend(additional_pins);
        self.base.requirements.extend(requirements);
        self.base.source = self
            .base
            .resolve_temperature(self.base.reference_temperature_k)?;
        self.static_recovery_properties = Some(properties);
        Ok(self)
    }

    /// Original independently selected property names, in retained branch order.
    pub const fn static_recovery_properties(&self) -> Option<&[String; N]> {
        self.static_recovery_properties.as_ref()
    }

    pub(super) fn static_recovery_law(
        &self,
        source: &ResolvedJ2Plasticity,
        law: &ChabocheJ2<N>,
        rate: Option<&ChabocheOverstressJ2<N>>,
        duration_s: Option<f64>,
    ) -> Result<Option<ChabocheStaticRecoveryJ2<N>>, FeedbackError> {
        let Some(properties) = &self.static_recovery_properties else {
            return Ok(None);
        };
        let duration = duration_s.ok_or(FeedbackError::Invalid(
            "source-defined static recovery requires the actual physical interval",
        ))?;
        let mut rates = [0.0; N];
        for (i, name) in properties.iter().enumerate() {
            rates[i] = source
                .material()
                .property(name)
                .ok_or(FeedbackError::Invalid(
                    "resolved source omitted an original static-recovery rate",
                ))?
                .value_si();
        }
        let prepared = match rate {
            Some(kinetic) => ChabocheStaticRecoveryJ2::with_overstress(kinetic.clone(), rates),
            None => ChabocheStaticRecoveryJ2::new(law.clone(), rates, duration),
        }
        .map_err(|error| FeedbackError::Source(error.into()))?;
        Ok(Some(prepared))
    }
}

impl<'a, const N: usize> CaloricPoint<'a, ThermalChabocheSource<'a, N>> {
    /// Select sourced dwell recovery before physical loading/time. The original
    /// uniform enthalpy owner then evolves temperature and native cyclic history
    /// together, using exactly one physical interval per accepted increment.
    pub fn with_static_recovery(
        mut self,
        selection: MaterialPropertySelection,
        properties: [String; N],
    ) -> Result<Self, FeedbackError> {
        if self.time_s != 0.0
            || self.strain != [0.0; 6]
            || self.stress_pa != [0.0; 6]
            || self.history.plastic.alpha != 0.0
            || self.history.plastic.plastic_strain != [0.0; 6]
            || self.history.backstress_pa != [[0.0; 6]; N]
        {
            return Err(FeedbackError::Invalid(
                "select static recovery before advancing the specimen",
            ));
        }
        self.binding = self.binding.with_static_recovery(selection, properties)?;
        self.source = self
            .binding
            .base()
            .resolve_temperature(self.state.temperature_k())?;
        Ok(self)
    }
}
