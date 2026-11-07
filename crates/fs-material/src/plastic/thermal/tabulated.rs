//! Original Cauchy hardening tables on the shared small-strain thermal owner.
//!
//! The native tabulated return, exact hardening storage, original caloric chart,
//! expansion and Prager history stay with their existing owners. Fixed strain
//! knots evolve affinely in temperature; their signed storage derivative enters
//! mechanical entropy and reversible heat. Yield heat and any selected Perzyna
//! loss remain physical dissipation, while integration loss is never heated.
//!
//! This family implements the existing thermal family boundary, so original
//! material-point feedback, spatial solids and layered sections retain the same
//! old history and physical interval throughout their coupled trials.

use super::expansion::ThermalHistoryContribution;
use super::family::{ThermoplasticFamilyError, ThermoplasticFamilyIncrement, ThermoplasticStep};
use super::*;
use crate::plastic::rate::{
    OVERSTRESS_EXPONENT_PROPERTY, OVERSTRESS_STRESS_SCALE_PROPERTY, OVERSTRESS_VISCOSITY_PROPERTY,
};
use crate::plastic::tabulated::rate::TabulatedOverstressJ2;
use crate::plastic::tabulated::source::CAUCHY_YIELD_CURVE_PROPERTY;
use crate::plastic::tabulated::thermal_source::{
    ResolvedThermalCauchyTabulatedHardening, ThermalHardeningCurveRequest,
    resolve_thermal_cauchy_tabulated_hardening,
};
use crate::plastic::tabulated::{TabulatedError, TabulatedJ2, TabulatedJ2Response};
use crate::state_point::YIELD_STRESS_PROPERTY;
use crate::{MaterialAdmissibility, Tangent6};
use fs_evidence::ModelCard;
use std::sync::Arc;

/// Same-card source family with the complete original Cauchy grid and receipts.
/// Cloning shares immutable table preparation and observations, never history.
#[derive(Clone, Debug)]
pub struct ThermalTabulatedSource<'a> {
    base: ThermoplasticSource<'a>,
    hardening: Arc<ResolvedThermalCauchyTabulatedHardening>,
}

/// One actual-temperature native law, with sourced thermal strain and optional
/// kinetics at exactly one physical duration. Repeated equilibrium trials must
/// reuse the same old J2State; accepted viscous history is inspected, not returned
/// again at the same endpoint to manufacture another interval of relaxation.
#[derive(Clone, Debug)]
pub struct ThermalTabulatedJ2 {
    source: ResolvedJ2Plasticity,
    law: TabulatedJ2,
    rate: Option<TabulatedOverstressJ2>,
    expansion_strain: f64,
}

impl<'a> ThermoplasticSource<'a> {
    /// Select the complete measured Cauchy yield radius on this same original
    /// card and caloric span. The separately supplied scalar initial yield must
    /// agree throughout its complete source line; energetic linear H must be an
    /// explicit zero. No scalar coefficient is inferred from a missing claim.
    /// Optional native linear/power kinetics should be selected before this call.
    pub fn with_tabulated_hardening(
        self,
        selection: ClaimSelection,
        maximum_source_samples: usize,
    ) -> Result<ThermalTabulatedSource<'a>, FeedbackError> {
        ThermalTabulatedSource::try_new(self, selection, maximum_source_samples)
    }
}

impl<'a> ThermalTabulatedSource<'a> {
    pub fn try_new(
        base: ThermoplasticSource<'a>,
        selection: ClaimSelection,
        maximum_source_samples: usize,
    ) -> Result<Self, FeedbackError> {
        if base.has_voce_hardening()
            || base.source.hardening_basis() == HardeningPropertyBasis::UniaxialTangentModulus
            || base.source.law().hardening != 0.0
            || base.hardening_temperature_slope_pa_k != 0.0
        {
            return Err(FeedbackError::Invalid(
                "thermal tables require explicit zero plastic-strain H and no competing Voce or total-strain tangent law",
            ));
        }
        let knots = base.caloric.curve().knots();
        let interval = [knots[0].temperature_k, knots[knots.len() - 1].temperature_k];
        let hardening = resolve_thermal_cauchy_tabulated_hardening(
            base.card,
            &base.context,
            &ThermalHardeningCurveRequest {
                temperature_axis: base.temperature_axis.clone(),
                reference_temperature_k: base.reference_temperature_k,
                temperature_range_k: interval,
                maximum_samples: maximum_source_samples,
                maximum_knots: maximum_source_samples,
            },
            selection,
        )
        .map_err(source_error)?;
        let property = base
            .source
            .material()
            .property(YIELD_STRESS_PROPERTY)
            .ok_or(FeedbackError::Invalid(
                "the original initial-yield claim is absent",
            ))?;
        let selected = property.answer().receipt.selected;
        let claim = base
            .card
            .claims()
            .claim(selected)
            .ok_or(FeedbackError::Invalid(
                "initial-yield claim left the shared material card",
            ))?;
        let lower = point_at(&base.context, &base.temperature_axis, interval[0])?;
        let upper = point_at(&base.context, &base.temperature_axis, interval[1])?;
        let line = base
            .card
            .claims()
            .query_line_typed(
                &claim.key,
                &base.temperature_axis,
                &lower,
                &upper,
                ClaimSelection::Pinned(selected),
            )
            .map_err(|error| FeedbackError::Support {
                property: claim.key.name().into(),
                reason: error.to_string(),
            })?;
        let samples: Vec<_> = line
            .knots
            .iter()
            .map(|(temperature, answer)| (*temperature, answer.evidence.value.value))
            .collect();
        let initial_yield = hardening.hardening().reference().initial_yield_stress_pa();
        let initial_slope = hardening.hardening().yield_stress_slopes_pa_k()[0];
        let scalar_slope = affine_temperature_slope(
            &samples,
            base.reference_temperature_k,
            base.source.law().yield_stress,
        )
        .map_err(FeedbackError::Invalid)?;
        if !agree(scalar_slope, initial_slope)
            || samples.iter().any(|&(temperature, value)| {
                !agree(
                    value,
                    initial_slope
                        .mul_add(temperature - base.reference_temperature_k, initial_yield),
                )
            })
            || !agree(base.source.law().yield_stress, initial_yield)
        {
            return Err(FeedbackError::Invalid(
                "scalar initial yield must agree with the original Cauchy table throughout the complete thermal source line",
            ));
        }
        finite(
            interval[1]
                * hardening
                    .hardening()
                    .maximum_absolute_storage_slope_j_m3_k(),
        )?;
        let result = Self {
            base,
            hardening: Arc::new(hardening),
        };
        // Every ordinate and segment difference is affine. The source already
        // admits both endpoint curves; the native law also admits return stiffness.
        for temperature in [
            interval[0],
            result.base.reference_temperature_k,
            interval[1],
        ] {
            result.native_at_temperature(temperature)?;
        }
        Ok(result)
    }

    pub const fn base(&self) -> &ThermoplasticSource<'a> {
        &self.base
    }

    /// Exact Cauchy observations, source decisions and affine storage preparation.
    pub fn hardening_source(&self) -> &ResolvedThermalCauchyTabulatedHardening {
        &self.hardening
    }

    fn native_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<(ResolvedJ2Plasticity, TabulatedJ2), FeedbackError> {
        let request = self.hardening.source().request();
        if !temperature_k.is_finite()
            || !(request.temperature_range_k[0]..=request.temperature_range_k[1])
                .contains(&temperature_k)
        {
            return Err(FeedbackError::Invalid(
                "thermal table temperature leaves its original admitted span",
            ));
        }
        let source = self.base.resolve_temperature(temperature_k)?;
        let table = self
            .hardening
            .hardening()
            .at_temperature_change(temperature_k - request.reference_temperature_k)
            .map_err(source_error)?;
        if !agree(source.law().yield_stress, table.initial_yield_stress_pa()) {
            return Err(FeedbackError::Invalid(
                "refreshed initial yield left the admitted thermal table",
            ));
        }
        let mut backbone = source.law().clone();
        // The explicitly selected table owns its entire absolute radius. The
        // original scalar receipt stays untouched; only roundoff differences
        // between its line evaluator and the admitted affine table are allowed.
        backbone.yield_stress = table.initial_yield_stress_pa();
        let law = TabulatedJ2::new(backbone, table)
            .map_err(|error| FeedbackError::Source(error.into()))?;
        Ok((source, law))
    }

    fn resolve_temperature(
        &self,
        temperature_k: f64,
        duration_s: Option<f64>,
    ) -> Result<ThermalTabulatedJ2, FeedbackError> {
        let (source, law) = self.native_at_temperature(temperature_k)?;
        let rate = if let Some(model) = self.base.overstress {
            let duration_s = duration_s.filter(|duration| positive(*duration)).ok_or(
                FeedbackError::Invalid(
                    "thermal tabulated overstress requires the actual positive physical duration",
                ),
            )?;
            let value = |name| {
                source
                    .material()
                    .property(name)
                    .map(|property| property.value_si())
                    .ok_or(FeedbackError::Invalid(
                        "resolved thermal table omitted an original kinetic coefficient",
                    ))
            };
            let viscosity = value(OVERSTRESS_VISCOSITY_PROPERTY)?;
            Some(
                match model {
                    OverstressSelection::Linear => {
                        TabulatedOverstressJ2::new_linear(law.clone(), viscosity, duration_s)
                    }
                    OverstressSelection::PowerLaw => TabulatedOverstressJ2::new_power_law(
                        law.clone(),
                        viscosity,
                        value(OVERSTRESS_STRESS_SCALE_PROPERTY)?,
                        value(OVERSTRESS_EXPONENT_PROPERTY)?,
                        duration_s,
                    ),
                }
                .map_err(|error| FeedbackError::Source(error.into()))?,
            )
        } else {
            None
        };
        let expansion_strain = self
            .base
            .expansion_strain_between(self.base.reference_temperature_k, temperature_k)?;
        Ok(ThermalTabulatedJ2 {
            source,
            law,
            rate,
            expansion_strain,
        })
    }

    /// Resolve an instantaneous law. Selected kinetics require the duration API.
    pub fn resolve_enthalpy(
        &self,
        h_j_kg: f64,
    ) -> Result<(EquilibriumPhaseState, ThermalTabulatedJ2), FeedbackError> {
        let phase = self
            .base
            .caloric
            .curve()
            .state_at_specific_enthalpy(h_j_kg)
            .map_err(FeedbackError::Caloric)?;
        let law = self.resolve_temperature(phase.temperature_k(), None)?;
        Ok((phase, law))
    }

    /// Resolve once for a physical interval; each Newton trial reuses old history.
    pub fn resolve_enthalpy_for_step(
        &self,
        h_j_kg: f64,
        duration_s: f64,
    ) -> Result<(EquilibriumPhaseState, ThermalTabulatedJ2), FeedbackError> {
        if !positive(duration_s) {
            return Err(FeedbackError::Invalid(
                "positive physical constitutive interval required",
            ));
        }
        let phase = self
            .base
            .caloric
            .curve()
            .state_at_specific_enthalpy(h_j_kg)
            .map_err(FeedbackError::Caloric)?;
        let law = self.resolve_temperature(phase.temperature_k(), Some(duration_s))?;
        Ok((phase, law))
    }

    fn validate_at_temperature(
        &self,
        temperature_k: f64,
        total_strain: &Voigt,
        state: &J2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FeedbackError> {
        if !positive(yield_tolerance_pa) {
            return Err(FeedbackError::Invalid(
                "positive finite yield allowance required",
            ));
        }
        let (_, law) = self.native_at_temperature(temperature_k)?;
        let mechanical = self.base.mechanical_strain(temperature_k, total_strain)?;
        if self.base.has_overstress() {
            law.backbone()
                .validate_overstress_history(&mechanical, state)
                .map_err(FeedbackError::Mechanical)?;
            law.hardening()
                .yield_stress_pa(state.alpha)
                .map_err(source_error)?;
            Ok(())
        } else {
            law.validate_committed_history(&mechanical, state, yield_tolerance_pa)
                .map_err(FeedbackError::Mechanical)
        }
    }
}

fn agree(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}

fn source_error(error: impl core::fmt::Display) -> FeedbackError {
    FeedbackError::Support {
        property: CAUCHY_YIELD_CURVE_PROPERTY.into(),
        reason: error.to_string(),
    }
}

impl ThermalTabulatedSource<'_> {
    /// Rate-independent native work and the original reciprocal caloric account.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_energy(
        &self,
        old_h_j_kg: f64,
        new_h_j_kg: f64,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<ThermoplasticIncrement, FeedbackError> {
        self.advance_for_interval(
            old_h_j_kg, new_h_j_kg, old_strain, old_state, strain, tolerances, None,
        )
    }

    /// Use the actual duration for the same native table and selected kinetics.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_energy_for_step(
        &self,
        old_h_j_kg: f64,
        new_h_j_kg: f64,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: f64,
    ) -> Result<ThermoplasticIncrement, FeedbackError> {
        if !positive(duration_s) {
            return Err(FeedbackError::Invalid(
                "positive physical constitutive interval required",
            ));
        }
        self.advance_for_interval(
            old_h_j_kg,
            new_h_j_kg,
            old_strain,
            old_state,
            strain,
            tolerances,
            Some(duration_s),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_for_interval(
        &self,
        old_h_j_kg: f64,
        new_h_j_kg: f64,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: Option<f64>,
    ) -> Result<ThermoplasticIncrement, FeedbackError> {
        let old = self
            .base
            .caloric
            .curve()
            .state_at_specific_enthalpy(old_h_j_kg)
            .map_err(FeedbackError::Caloric)?;
        let new = self
            .base
            .caloric
            .curve()
            .state_at_specific_enthalpy(new_h_j_kg)
            .map_err(FeedbackError::Caloric)?;
        Ok(self
            .increment_at_temperature(
                old.temperature_k(),
                new.temperature_k(),
                old_strain,
                old_state,
                strain,
                tolerances,
                duration_s,
            )?
            .1)
    }

    #[allow(clippy::too_many_arguments)]
    fn increment_at_temperature(
        &self,
        old_temperature_k: f64,
        temperature_k: f64,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: Option<f64>,
    ) -> Result<(ResolvedJ2Plasticity, ThermoplasticIncrement), FeedbackError> {
        self.validate_at_temperature(
            old_temperature_k,
            old_strain,
            old_state,
            tolerances.yield_pa,
        )?;
        let endpoint = self.resolve_temperature(temperature_k, duration_s)?;
        let old_mechanical = self.base.mechanical_strain(old_temperature_k, old_strain)?;
        let mechanical_strain = endpoint.mechanical_strain(strain);
        let mechanical = match &endpoint.rate {
            Some(rate) => {
                rate.advance_with_energy(&old_mechanical, old_state, &mechanical_strain, tolerances)
            }
            None => endpoint.law.advance_with_energy_from_admitted_history(
                &old_mechanical,
                old_state,
                &mechanical_strain,
                tolerances,
            ),
        }
        .map_err(FeedbackError::Mechanical)?;
        let mut contribution = ThermalHistoryContribution::default();
        let hardening = self.hardening.hardening();
        if hardening.maximum_absolute_storage_slope_j_m3_k() != 0.0 {
            let old_derivative = hardening
                .hardening_energy_temperature_derivative_j_m3_k(old_state.alpha)
                .map_err(source_error)?;
            let new_derivative = hardening
                .hardening_energy_temperature_derivative_j_m3_k(mechanical.state.alpha)
                .map_err(source_error)?;
            contribution.variable_storage = true;
            contribution.old_entropy_j_m3_k = -old_derivative;
            contribution.new_entropy_j_m3_k = -new_derivative;
            contribution.coefficient_storage_change_j_m3 =
                product(temperature_k - old_temperature_k, old_derivative)?;
        }
        let increment = self.base.finish_increment_at_temperature(
            old_temperature_k,
            temperature_k,
            old_strain,
            old_state,
            strain,
            &endpoint.source,
            mechanical,
            contribution,
            tolerances,
        )?;
        Ok((endpoint.source, increment))
    }
}

impl ThermalTabulatedJ2 {
    /// Untouched original mechanical, kinetic and caloric source bundle.
    pub const fn resolved(&self) -> &ResolvedJ2Plasticity {
        &self.source
    }

    /// Native bounded table at this actual temperature, before optional kinetics.
    pub const fn tabulated_law(&self) -> &TabulatedJ2 {
        &self.law
    }

    pub const fn rate_law(&self) -> Option<&TabulatedOverstressJ2> {
        self.rate.as_ref()
    }

    fn mechanical_strain(&self, strain: &Voigt) -> Voigt {
        core::array::from_fn(|i| strain[i] - if i < 3 { self.expansion_strain } else { 0.0 })
    }

    /// Checked original tensor return and consistent tangent at fixed old history.
    pub fn response(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<TabulatedJ2Response, TabulatedError> {
        let mechanical = self.mechanical_strain(strain);
        match &self.rate {
            Some(rate) => rate.response(&mechanical, state),
            None => self.law.response(&mechanical, state),
        }
    }
}

impl SmallStrainLaw for ThermalTabulatedJ2 {
    type State = J2State;

    fn initial_state(&self) -> J2State {
        self.law.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &J2State) -> Voigt {
        let strain = self.mechanical_strain(strain);
        match &self.rate {
            Some(rate) => rate.stress(&strain, state),
            None => self.law.stress(&strain, state),
        }
    }
    fn tangent(&self, strain: &Voigt, state: &J2State) -> Tangent6 {
        let strain = self.mechanical_strain(strain);
        match &self.rate {
            Some(rate) => rate.tangent(&strain, state),
            None => self.law.tangent(&strain, state),
        }
    }
    fn update_state(&self, strain: &Voigt, state: &J2State) -> J2State {
        let strain = self.mechanical_strain(strain);
        match &self.rate {
            Some(rate) => rate.update_state(&strain, state),
            None => self.law.update_state(&strain, state),
        }
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        match &self.rate {
            Some(rate) => rate.admissibility(),
            None => self.law.admissibility(),
        }
    }
    fn card(&self) -> ModelCard {
        match &self.rate {
            Some(rate) => rate.card(),
            None => self.law.card(),
        }
    }
}

impl<'a> ThermoplasticSourceFamily<'a> for ThermalTabulatedSource<'a> {
    type Material = ThermalTabulatedJ2;

    fn thermal_source(&self) -> &ThermoplasticSource<'a> {
        &self.base
    }

    fn resolve_step(
        &self,
        old_h: f64,
        new_h: f64,
        duration_s: f64,
    ) -> Result<ThermoplasticStep<Self::Material>, ThermoplasticFamilyError> {
        self.base
            .caloric
            .curve()
            .state_at_specific_enthalpy(old_h)
            .map_err(FeedbackError::Caloric)?;
        let (phase, law) = self.resolve_enthalpy_for_step(new_h, duration_s)?;
        Ok(ThermoplasticStep {
            temperature_k: phase.temperature_k(),
            material: law.resolved().clone(),
            law,
        })
    }

    fn checked_increment(
        &self,
        old_h: f64,
        new_h: f64,
        old_strain: &Voigt,
        old: &J2State,
        new_strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: f64,
    ) -> Result<ThermoplasticFamilyIncrement<J2State>, ThermoplasticFamilyError> {
        let thermal = self.advance_with_energy_for_step(
            old_h, new_h, old_strain, old, new_strain, tolerances, duration_s,
        )?;
        Ok(ThermoplasticFamilyIncrement {
            state: thermal.mechanical.state.clone(),
            thermal,
        })
    }

    fn checked_temperature_increment(
        &self,
        old_temperature_k: f64,
        temperature_k: f64,
        old_strain: &Voigt,
        old: &J2State,
        new_strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: f64,
    ) -> Result<
        (ResolvedJ2Plasticity, ThermoplasticFamilyIncrement<J2State>),
        ThermoplasticFamilyError,
    > {
        if !positive(duration_s) {
            return Err(
                FeedbackError::Invalid("positive physical constitutive interval required").into(),
            );
        }
        let (source, thermal) = self.increment_at_temperature(
            old_temperature_k,
            temperature_k,
            old_strain,
            old,
            new_strain,
            tolerances,
            Some(duration_s),
        )?;
        Ok((
            source,
            ThermoplasticFamilyIncrement {
                state: thermal.mechanical.state.clone(),
                thermal,
            },
        ))
    }

    fn validate_accepted(
        &self,
        h: f64,
        strain: &Voigt,
        state: &J2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), ThermoplasticFamilyError> {
        let phase = self
            .base
            .caloric
            .curve()
            .state_at_specific_enthalpy(h)
            .map_err(FeedbackError::Caloric)?;
        self.validate_at_temperature(phase.temperature_k(), strain, state, yield_tolerance_pa)
            .map_err(Into::into)
    }
}

impl<'a> CaloricPoint<'a> {
    /// Attach original table evidence before loading or physical time advances.
    /// Select optional source-defined kinetics first; the existing feedback loop
    /// supplies its actual duration and receives physical/reversible heat once.
    pub fn with_tabulated_hardening(
        self,
        selection: ClaimSelection,
        maximum_source_samples: usize,
    ) -> Result<CaloricPoint<'a, ThermalTabulatedSource<'a>>, FeedbackError> {
        if self.time_s != 0.0
            || self.strain != [0.0; 6]
            || self.stress_pa != [0.0; 6]
            || self.history.alpha != 0.0
            || self.history.plastic_strain != [0.0; 6]
        {
            return Err(FeedbackError::Invalid(
                "select thermal tabulated hardening before advancing the specimen",
            ));
        }
        let binding = self
            .binding
            .with_tabulated_hardening(selection, maximum_source_samples)?;
        let source = binding.base.reference().clone();
        Ok(CaloricPoint {
            binding,
            source,
            strain: self.strain,
            stress_pa: self.stress_pa,
            history: self.history,
            state: self.state,
            mass_kg: self.mass_kg,
            volume_m3: self.volume_m3,
            surface_area_m2: self.surface_area_m2,
            time_s: self.time_s,
            source_lifetime: core::marker::PhantomData,
        })
    }
}

impl<'a> CaloricPoint<'a, ThermalTabulatedSource<'a>> {
    pub fn tabulated_hardening_source(&self) -> &ResolvedThermalCauchyTabulatedHardening {
        self.binding.hardening_source()
    }
}
