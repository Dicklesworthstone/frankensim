//! Reversible heat on the original uniform string/transport transaction.
//!
//! Mechanics advances once at the accepted temperature. Only the endpoint heat
//! exchange is iterated, using the existing bounded scalar feedback owner.

use super::{
    AcousticRealizeError, LinearMaterialStringRuntime, MaterialStringFrame,
    ThermalMaterialStringRuntime, ThermalStringError, string_mode_omega,
};
use crate::cavity_phs::piston::thermal::feedback::{self, ScalarError};
use fs_material::phase::{
    EquilibriumPhaseState, SolidLiquidPhase, UniformEnthalpyStep, UniformEnthalpyStepInput,
};

#[derive(Default)]
pub(super) struct ThermalStringExchange {
    pub axial_free_energy_change_j: f64,
    pub mechanical_entropy_energy_change_j: f64,
    pub modal_parameter_work_j: f64,
    pub reversible_heat_j: f64,
    pub roundoff_tolerance_j: f64,
}

fn check(cx: &fs_exec::Cx<'_>) -> Result<(), ThermalStringError> {
    cx.checkpoint()
        .map_err(|_| ThermalStringError::Acoustic(AcousticRealizeError::Cancelled))
}

impl ThermalMaterialStringRuntime {
    /// Accepted mechanical internal energy [J], excluding the independent
    /// caloric chart. Includes kinetic/modal energy and, when thermoelastic,
    /// fixed-support axial free energy plus `T*s` from both elastic stores.
    /// This is the storage to use in a combined electrical/solid/gas balance;
    /// modal energy alone omits the reversible temperature contribution.
    pub fn mechanical_internal_energy_j(&self) -> Result<f64, ThermalStringError> {
        let modal: f64 = self.mechanical.modes().iter().zip(self.mechanical.states())
            .map(|(mode, state)| 0.5 * (state.velocity_m_sqrt_kg_per_s.powi(2)
                + (mode.angular_frequency_rad_s * state.displacement_m_sqrt_kg).powi(2)))
            .sum();
        let energy = if let Some(law) = &self.thermoelastic {
            let state = law.energy(&self.mechanical, self.thermal.temperature_k())?;
            modal + state.axial_free_energy_j + self.thermal.temperature_k() * state.entropy_j_k
        } else { modal };
        if !energy.is_finite() {
            return Err(ThermalStringError::Admission("mechanical internal energy is not finite"));
        }
        Ok(energy)
    }

    /// Every proposal uses the original thermal state and the SAME completed
    /// mechanical substep. Source errors keep their concrete enclosing type.
    pub(super) fn propose_thermal_step<R, E: From<ThermalStringError>>(
        &self,
        cx: &fs_exec::Cx<'_>,
        mechanical: &LinearMaterialStringRuntime,
        vibration: &MaterialStringFrame,
        mut transport: impl FnMut(UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<R>, E>,
    ) -> Result<UniformEnthalpyStep<R>, E> {
        check(cx)?;
        let specimen = self.mechanical.specimen();
        let radius = specimen.radius_m();
        let length = specimen.string().length_m;
        let base = UniformEnthalpyStepInput {
            curve: &self.curve,
            initial: self.thermal,
            mass_kg: specimen.mass_kg(),
            volume_m3: specimen.area_m2() * length,
            surface_area_m2: 2.0 * core::f64::consts::PI * radius * (length + radius),
            internal_heat_j: vibration.dissipation.material_heat_j,
            duration_s: mechanical.member.model.sample_period_s(),
        };
        let mut propose = |reversible_heat: f64| {
            check(cx)?;
            let internal_heat_j = base.internal_heat_j + reversible_heat;
            if !internal_heat_j.is_finite() {
                return Err(E::from(ThermalStringError::Admission(
                    "material dissipation and reversible heat overflow",
                )));
            }
            let proposal = transport(UniformEnthalpyStepInput { internal_heat_j, ..base })?;
            self.check_thermal_proposal(internal_heat_j, &proposal)?;
            check(cx)?;
            Ok(proposal)
        };
        let Some(law) = &self.thermoelastic else {
            return propose(0.0);
        };
        // Modal motion changes the fixed-coordinate heat capacity as well as
        // entropy. Recheck the actual proposed vibration before any heat solve.
        law.admit_capacity(cx, mechanical, &self.curve)?;
        let controls = law.controls();
        let mut proposal = feedback::solve_scalar(
            controls.reversible_heat_bracket_j,
            controls.energy_tolerance_j,
            controls.maximum_iterations,
            |applied, _, _| {
                let proposed = propose(applied)?;
                let (_, exchange) = self.prepare_thermal_rebind(cx, mechanical, proposed.state)?;
                let residual = applied - exchange.reversible_heat_j;
                Ok::<_, E>((proposed, residual))
            },
            &mut || cx.checkpoint().is_err(),
        ).map_err(|error| match error {
            ScalarError::Evaluation(error) => error,
            ScalarError::Invalid(what) => E::from(ThermalStringError::Admission(what)),
            ScalarError::Cancelled => E::from(ThermalStringError::Acoustic(AcousticRealizeError::Cancelled)),
            ScalarError::OutsideBracket { lower, upper } => E::from(ThermalStringError::ThermoelasticBracket {
                lower_residual_j: lower, upper_residual_j: upper,
            }),
            ScalarError::NotConverged { iterations, residual } => E::from(ThermalStringError::ThermoelasticConvergence {
                iterations, residual_j: residual,
            }),
        })?;
        // Keep the transport's independent tolerance AND the heat-source solve
        // allowance. Neither one is a bound on the first-order partition error.
        proposal.energy_residual_tolerance_j += controls.energy_tolerance_j;
        if !proposal.energy_residual_tolerance_j.is_finite() {
            return Err(ThermalStringError::Admission("combined thermal solve allowance overflows").into());
        }
        Ok(proposal)
    }

    fn check_thermal_proposal<R>(
        &self,
        internal_heat_j: f64,
        proposal: &UniformEnthalpyStep<R>,
    ) -> Result<(), ThermalStringError> {
        let checked = self.curve
            .state_at_specific_enthalpy(proposal.state.specific_enthalpy_j_kg())
            .map_err(ThermalStringError::Phase)?;
        if checked != proposal.state || !proposal.external_heat_j.is_finite()
            || !proposal.energy_residual_tolerance_j.is_finite()
            || proposal.energy_residual_tolerance_j < 0.0 {
            return Err(ThermalStringError::Admission(
                "thermal transport returned a foreign state, invalid heat or invalid residual budget",
            ));
        }
        let mass = self.mechanical.specimen().mass_kg();
        let old_h = self.thermal.specific_enthalpy_j_kg();
        let new_h = proposal.state.specific_enthalpy_j_kg();
        let change = mass * (new_h - old_h);
        let residual = change - internal_heat_j - proposal.external_heat_j;
        let roundoff = 32.0 * f64::EPSILON * (mass * (old_h.abs() + new_h.abs())
            + internal_heat_j.abs() + proposal.external_heat_j.abs());
        if !residual.is_finite() || !roundoff.is_finite()
            || residual.abs() > roundoff + proposal.energy_residual_tolerance_j {
            return Err(ThermalStringError::Admission(
                "thermal transport violates its own caloric heat balance",
            ));
        }
        Ok(())
    }

    /// Rebind through the original acoustic owner, using only its calculated
    /// modal work as the admission budget. No externally invented work or
    /// rescaling of modal coordinates is allowed.
    pub(super) fn prepare_thermal_rebind(
        &self,
        cx: &fs_exec::Cx<'_>,
        mechanical: &LinearMaterialStringRuntime,
        thermal: EquilibriumPhaseState,
    ) -> Result<(LinearMaterialStringRuntime, ThermalStringExchange), ThermalStringError> {
        check(cx)?;
        if thermal.phase_curve_identity() != self.curve.identity()
            || thermal.phase() != SolidLiquidPhase::Solid
            || (thermal.temperature_k() - self.thermal.temperature_k()).abs()
                > self.maximum_temperature_increment_k {
            return Err(ThermalStringError::Admission(
                "thermal step left its solid chart or exceeded its temperature-increment budget",
            ));
        }
        let specimen = self.specimen_at(thermal)?;
        let new_string = specimen.string();
        // Exactly the same native frequency and difference-of-squares used
        // by rebind; a zero-work ordinary string keeps its original strict gate.
        let work_budget = if self.thermoelastic.is_some() {
            mechanical.modes().iter().zip(mechanical.states()).enumerate()
                .map(|(index, (mode, state))| {
                    let old = mode.angular_frequency_rad_s * state.displacement_m_sqrt_kg;
                    let new = string_mode_omega(&new_string, index + 1) * state.displacement_m_sqrt_kg;
                    0.5 * (new - old) * (new + old)
                }).sum::<f64>().abs()
        } else { 0.0 };
        let mut rebound = mechanical.clone();
        // Share the sample's single new epoch, including the final representable one.
        rebound.epoch = self.mechanical.epoch();
        let update = rebound.rebind(cx, specimen, work_budget)
            .map_err(ThermalStringError::Acoustic)?;
        let Some(law) = &self.thermoelastic else {
            return Ok((rebound, ThermalStringExchange::default()));
        };
        let before = law.energy(&self.mechanical, self.thermal.temperature_k())?;
        let after = law.energy(&rebound, thermal.temperature_k())?;
        let before_ts = self.thermal.temperature_k() * before.entropy_j_k;
        let after_ts = thermal.temperature_k() * after.entropy_j_k;
        let axial_free_energy_change_j = after.axial_free_energy_j - before.axial_free_energy_j;
        let mechanical_entropy_energy_change_j = after_ts - before_ts;
        let modal_parameter_work_j = update.parameter_work_j;
        let reversible_heat_j = -axial_free_energy_change_j - modal_parameter_work_j
            - mechanical_entropy_energy_change_j;
        let scale = before.axial_free_energy_j.abs() + after.axial_free_energy_j.abs()
            + before.modal_free_energy_j.abs() + after.modal_free_energy_j.abs()
            + before_ts.abs() + after_ts.abs() + modal_parameter_work_j.abs();
        let roundoff_tolerance_j = 128.0 * f64::EPSILON * scale;
        if [axial_free_energy_change_j, mechanical_entropy_energy_change_j,
            modal_parameter_work_j, reversible_heat_j, roundoff_tolerance_j]
            .iter().any(|value| !value.is_finite()) {
            return Err(ThermalStringError::Admission("thermoelastic storage or heat is not finite"));
        }
        check(cx)?;
        Ok((rebound, ThermalStringExchange {
            axial_free_energy_change_j, mechanical_entropy_energy_change_j,
            modal_parameter_work_j, reversible_heat_j, roundoff_tolerance_j,
        }))
    }
}
