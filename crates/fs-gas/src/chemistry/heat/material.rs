//! A finite caloric body exchanges heat inside the original chemical stages.
//!
//! The body follows its original equilibrium h chart at an explicitly held
//! pressure. Its mass is fixed. Heat changes H; the original material owner
//! reports expansion work and the corresponding U change. Gas and body heat
//! cancel once, including through sharp latent plateaus. This is a homogeneous
//! body approximation, not a spatial wall, moving boundary or contact model.

use super::*;
use crate::chemistry::isobaric::{IsobaricKineticsView, IsobaricReactor};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState,
    isobaric::IsobaricHeatStep,
};

/// Fixed-mass body with its exact original caloric source and current state.
///
/// The caller supplies the pressure at which the chart applies and a total
/// conductance G [W/K]. Density changes never change mass or gas volume. Internal
/// gradients must be negligible for this homogeneous approximation; a smaller
/// timestep cannot justify a spatially invalid lumped model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteHeatBody<'m> {
    curve: &'m EquilibriumEnthalpyPhaseCurve,
    state: EquilibriumPhaseState,
    mass: Mass,
    pressure: Pressure,
    conductance: ThermalConductance,
}

impl<'m> FiniteHeatBody<'m> {
    /// Bound the source traversal used by each exact thermal entropy integral.
    pub const MAXIMUM_KNOTS: usize = 4096;

    pub fn new(
        curve: &'m EquilibriumEnthalpyPhaseCurve,
        mass: Mass,
        pressure: Pressure,
        specific_enthalpy_j_kg: f64,
        conductance: ThermalConductance,
    ) -> Result<Self, IsochoricKineticsError> {
        if !positive(mass.value())
            || !pressure.value().is_finite()
            || pressure.value() < 0.0
            || !conductance.value().is_finite()
            || conductance.value() < 0.0
            || curve.knots().len() > Self::MAXIMUM_KNOTS
        {
            return Err(IsochoricKineticsError::Invalid(
                "invalid finite thermal body or source knot budget",
            ));
        }
        let state = curve.state_at_specific_enthalpy(specific_enthalpy_j_kg)?;
        // Admit actual mass, volume and pressure-work representation before time.
        curve.advance_isobaric_heat(state, mass.value(), pressure.value(), 0.0)?;
        Ok(Self {
            curve,
            state,
            mass,
            pressure,
            conductance,
        })
    }

    pub const fn curve(self) -> &'m EquilibriumEnthalpyPhaseCurve {
        self.curve
    }
    pub const fn state(self) -> EquilibriumPhaseState {
        self.state
    }
    pub const fn mass(self) -> Mass {
        self.mass
    }
    pub const fn pressure(self) -> Pressure {
        self.pressure
    }
    pub const fn conductance(self) -> ThermalConductance {
        self.conductance
    }

    pub(super) fn sample(
        self,
        gas_heat_j: f64,
    ) -> Result<EnthalpyPhaseKnot, IsochoricKineticsError> {
        let increment = -gas_heat_j / self.mass.value();
        if !increment.is_finite() || (gas_heat_j != 0.0 && increment == 0.0) {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable finite-body specific heat",
            ));
        }
        Ok(self
            .curve
            .values_at_specific_enthalpy(self.state.specific_enthalpy_j_kg() + increment)?)
    }

    fn heat_report(self, gas_heat_j: f64) -> Result<IsobaricHeatStep, IsochoricKineticsError> {
        Ok(self.curve.advance_isobaric_heat(
            self.state,
            self.mass.value(),
            self.pressure.value(),
            -gas_heat_j,
        )?)
    }

    pub(super) fn admit_heat(
        self,
        gas_heat_j: f64,
        tolerance: f64,
    ) -> Result<(), IsochoricKineticsError> {
        if self.heat_report(gas_heat_j)?.heat_residual_j.abs() > tolerance {
            return Err(IsochoricKineticsError::Invalid(
                "finite-body enthalpy datum cannot resolve the energy tolerance",
            ));
        }
        Ok(())
    }

    pub(super) fn entropy_transfer(
        self,
        initial: Sum,
        final_heat: Sum,
        volume: f64,
    ) -> Result<f64, IsochoricKineticsError> {
        let initial = self.sample(initial.value() * volume)?;
        let final_state = self.sample(final_heat.value() * volume)?;
        Ok(-self.mass.value()
            * self.curve.specific_entropy_change_j_kg_k(
                initial.specific_enthalpy_j_kg,
                final_state.specific_enthalpy_j_kg,
            )?)
    }

    pub(super) fn additional_error(
        self,
        base: Sum,
        full: Sum,
        fine: Sum,
        volume: f64,
        controls: IsochoricKineticsSolve,
    ) -> Result<f64, IsochoricKineticsError> {
        let base = self.sample(base.value() * volume)?;
        let full = self.sample(full.value() * volume)?;
        let fine = self.sample(fine.value() * volume)?;
        let temperature = (fine.temperature_k - full.temperature_k).abs()
            / (controls.relative_tolerance * fine.temperature_k.max(base.temperature_k));
        // Fraction is dimensionless on [0,1]; it still changes on a T plateau.
        let phase = (fine.liquid_mass_fraction - full.liquid_mass_fraction).abs()
            / controls.relative_tolerance;
        Ok(temperature.max(phase))
    }

    // q = fixed + a*(Tb(h)-Tg), h = h0 - V*q/m. On each original segment,
    // h + beta*Tb(h) is affine and strictly increasing, including a plateau.
    // Locate that segment by binary search, then eliminate h analytically.
    // No added Newton coordinate, nested iterative root or apparent Cp exists.
    pub(super) fn stage_heat(
        self,
        fixed: Sum,
        volume: f64,
        capacity: f64,
        temperature: f64,
    ) -> Result<(Sum, f64), IsochoricKineticsError> {
        if capacity == 0.0 {
            self.sample(fixed.value() * volume)?;
            return Ok((fixed, 0.0));
        }
        let beta = capacity * volume / self.mass.value();
        let forcing = fixed.value() * volume / self.mass.value();
        if !beta.is_finite() || !forcing.is_finite() {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable finite-body thermal stage",
            ));
        }
        let residual = |knot: &EnthalpyPhaseKnot| {
            let mut value = Sum::default();
            value.add(knot.specific_enthalpy_j_kg - self.state.specific_enthalpy_j_kg());
            value.add(forcing);
            value.add(beta * (knot.temperature_k - temperature));
            value.value()
        };
        let knots = self.curve.knots();
        let upper = knots
            .partition_point(|knot| residual(knot) <= 0.0)
            .clamp(1, knots.len() - 1);
        let lower = knots[upper - 1];
        let slope = (knots[upper].temperature_k - lower.temperature_k)
            / (knots[upper].specific_enthalpy_j_kg - lower.specific_enthalpy_j_kg);
        let denominator = 1.0 + beta * slope;
        if !slope.is_finite() || !positive(denominator) {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable finite-body caloric slope or implicit denominator",
            ));
        }
        let h = lower.specific_enthalpy_j_kg - residual(&lower) / denominator;
        let body_temperature = self.curve.temperature_at_specific_enthalpy(h)?;
        let mut heat = fixed;
        heat.add(capacity * (body_temperature - temperature));
        let feedback = capacity / denominator;
        if !heat.value().is_finite() || !feedback.is_finite() || feedback < 0.0 {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable finite-body thermal feedback",
            ));
        }
        Ok((heat, feedback))
    }

    pub(super) fn stage_temperature_interval(
        self,
        fixed: Sum,
        volume: f64,
        capacity: f64,
        minimum: f64,
        maximum: f64,
    ) -> Result<(f64, f64), IsochoricKineticsError> {
        if capacity == 0.0 {
            self.sample(fixed.value() * volume)?;
            return Ok((minimum, maximum));
        }
        let at = |knot: &EnthalpyPhaseKnot| {
            let mut heat = fixed;
            heat.add(
                (knot.specific_enthalpy_j_kg - self.state.specific_enthalpy_j_kg())
                    * (self.mass.value() / volume),
            );
            knot.temperature_k + heat.value() / capacity
        };
        let knots = self.curve.knots();
        // These are gas temperatures whose implicit body state is exactly at
        // its source bounds. Gas caloric bracketing must not probe beyond them.
        let lower = minimum.max(at(&knots[0]));
        let upper = maximum.min(at(&knots[knots.len() - 1]));
        if !lower.is_finite() || !upper.is_finite() || lower > upper {
            return Err(IsochoricKineticsError::Invalid(
                "gas and finite-body caloric domains do not overlap",
            ));
        }
        Ok((lower, upper))
    }

    fn finish<G>(
        self,
        gas: G,
        heat_j: f64,
        tolerance: f64,
    ) -> Result<FiniteHeatBodyStep<'m, G>, IsochoricKineticsError> {
        self.admit_heat(heat_j, tolerance)?;
        let report = self.heat_report(heat_j)?;
        let entropy_change_j_k = self.mass.value()
            * self.curve.specific_entropy_change_j_kg_k(
                self.state.specific_enthalpy_j_kg(),
                report.state.specific_enthalpy_j_kg(),
            )?;
        Ok(FiniteHeatBodyStep {
            gas,
            body: Self {
                state: report.state,
                ..self
            },
            report,
            entropy_change_j_k,
        })
    }
}

/// Both complete endpoints publish together. `body()` can continue the original
/// source at the next interval; no body state changes on a failed gas advance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteHeatBodyStep<'m, G> {
    gas: G,
    body: FiniteHeatBody<'m>,
    report: IsobaricHeatStep,
    entropy_change_j_k: f64,
}
impl<'m, G> FiniteHeatBodyStep<'m, G> {
    pub const fn gas(&self) -> G
    where
        G: Copy,
    {
        self.gas
    }
    pub const fn body(&self) -> FiniteHeatBody<'m> {
        self.body
    }
    /// Actual material H/U changes, density, volumes, p*dV and heat residual.
    pub const fn body_report(&self) -> &IsobaricHeatStep {
        &self.report
    }
    pub const fn body_entropy_change_j_k(&self) -> f64 {
        self.entropy_change_j_k
    }
    pub fn into_parts(self) -> (G, FiniteHeatBody<'m>, IsobaricHeatStep) {
        (self.gas, self.body, self.report)
    }
}

impl IsochoricReactor {
    /// Couple a fixed-volume reacting parcel and the original finite body.
    /// Ugas + Hbody closes before body expansion work is removed; equivalently
    /// delta Ugas + delta Ubody + body p*dV = 0. Reaction energy enters U once.
    /// Both BE/SDIRK2 include body cooling/melting in every chemical stage.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_finite_body<'a, 'm>(
        &'a mut self,
        amounts: &[Amount],
        volume: Volume,
        initial_temperature: Temperature,
        duration: Time,
        body: FiniteHeatBody<'m>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteHeatBodyStep<'m, IsochoricKineticsView<'a>>, IsochoricKineticsError> {
        let tolerance = self.controls().energy_tolerance_j;
        let gas = advance_prepared_constraint(
            self,
            amounts,
            None,
            volume,
            None,
            initial_temperature,
            duration,
            Some(HeatExchange::Body(body)),
            &mut cancelled,
        )?;
        let heat = gas.balance().heat_input_j;
        let result = body.finish(gas, heat, tolerance)?;
        poll(&mut cancelled)?;
        Ok(result)
    }
}

impl IsobaricReactor {
    /// Couple a held-pressure reacting parcel and a finite original caloric body.
    /// Hgas + Hbody closes. Each original owner reports its own outward p*dV;
    /// those two work terms close the enclosing internal-energy first law.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_finite_body<'a, 'm>(
        &'a mut self,
        amounts: &[Amount],
        pressure: Pressure,
        initial_temperature: Temperature,
        duration: Time,
        body: FiniteHeatBody<'m>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteHeatBodyStep<'m, IsobaricKineticsView<'a>>, IsochoricKineticsError> {
        let tolerance = self.controls().energy_tolerance_j;
        let gas = self.advance(
            amounts,
            None,
            pressure,
            initial_temperature,
            duration,
            Some(HeatExchange::Body(body)),
            &mut cancelled,
        )?;
        let heat = gas.balance().heat_input_j;
        let result = body.finish(gas, heat, tolerance)?;
        poll(&mut cancelled)?;
        Ok(result)
    }
}
