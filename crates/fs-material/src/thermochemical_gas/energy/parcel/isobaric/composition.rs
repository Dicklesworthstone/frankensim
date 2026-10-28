//! Composition and heat at one constant-pressure material endpoint.
//!
//! Total formation enthalpy is part of H, not an extra heat source. The
//! original root solver reconstructs h_new(T,p)-h_old(T0,p)=Q/m; the same
//! endpoint EOS supplies V1. Independent Delta U and actual p(V1-V0) close
//! the first law. Stoichiometry and kinetic evolution belong to the caller.
//! Reference: Cantera 3.2, constant-pressure-reactor energy equation.

use super::ClosedGasParcel;
use super::super::composition::GasCompositionStep;
use crate::thermochemical_gas::{ThermochemicalGasModel, failure};
use crate::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, map_root, poll, positive, root, updated_point,
};
use fs_qty::{Pressure, Temperature};

/// Source-resolved composition update and its independent enthalpy account.
#[derive(Clone, Debug, PartialEq)]
pub struct GasIsobaricCompositionStep {
    /// Original gas and constant external load pressure [Pa].
    pub pressure_pa: f64,
    /// Original extensive internal-energy, heat and actual volume-work account.
    pub energy: GasCompositionStep,
    /// Actual m(h_new-h_old), including chemical formation references [J].
    pub enthalpy_change_j: f64,
    /// Actual enthalpy change minus supplied heat [J].
    pub enthalpy_residual_j: f64,
}

impl ClosedGasParcel {
    /// Reconstruct a new composition at the CURRENT pressure and fixed mass.
    ///
    /// The caller owns species/element/charge conservation and reaction rates.
    /// The two models must share a standard-state reference convention. No
    /// additional reaction heat, constant Cp, temperature average or frozen
    /// density is introduced. Expansion/contraction work is p*(V_new-V_old),
    /// not an independently adjustable input. Shaft work and piston inertia
    /// are outside this quasistatic pressure constraint.
    ///
    /// Only the NEW continuous caloric branch must cover the solve bracket.
    /// There is no fictitious intermediate at new composition and old energy
    /// or temperature. Full source/transport data resolve once at the endpoint.
    /// Identical models delegate to the reference-safe frozen-composition heat
    /// owner. Changed compositions retain absolute formation references; if
    /// their differences cannot close the specified absolute balances, refuse
    /// rather than replacing chemical energy by sensible heat.
    ///
    /// Source, convention, work, budget, arithmetic or cancellation refusal
    /// preserves the complete accepted parcel. No entropy-production assertion
    /// is made here: the kinetic/thermal owner must audit its physical process.
    pub fn replace_model_with_heat_at_constant_pressure(
        &mut self,
        model: ThermochemicalGasModel,
        heat_j: f64,
        controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasIsobaricCompositionStep, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !heat_j.is_finite() {
            return Err(GasEnergyError::Invalid("finite isobaric composition heat required"));
        }
        let initial = self.parcel.state.clone();
        if model == self.parcel.model {
            let step = self.add_heat_at_constant_pressure(heat_j, controls, &mut cancelled)?;
            return Ok(GasIsobaricCompositionStep {
                pressure_pa: step.pressure_pa,
                energy: GasCompositionStep {
                    initial, final_state: self.parcel.state.clone(), mass_kg: step.mass_kg,
                    initial_volume_m3: step.initial_volume_m3, volume_m3: step.final_volume_m3,
                    supplied_heat_j: heat_j, work_by_gas_j: step.work_by_gas_j,
                    stored_energy_change_j: step.stored_energy_change_j,
                    energy_residual_j: step.energy_residual_j, iterations: step.iterations,
                },
                enthalpy_change_j: step.enthalpy_change_j,
                enthalpy_residual_j: step.enthalpy_residual_j,
            });
        }
        let old = &self.parcel;
        let convention = old.model.thermodynamics().components()[0].convention();
        if model.thermodynamics().components().iter().any(|s| s.convention() != convention) {
            return Err(GasEnergyError::Invalid(
                "isobaric composition changes require one common thermochemical reference convention",
            ));
        }
        let (mass, volume0, pressure) = (old.mass_kg, old.volume_m3, initial.state().pressure);
        let q = heat_j / mass;
        let tolerance_j = mass * controls.energy_tolerance_j_kg;
        if !positive(mass) || !q.is_finite() || (heat_j != 0.0 && q == 0.0)
            || !positive(tolerance_j)
        {
            return Err(GasEnergyError::Invalid("isobaric composition energy scale is unrepresentable"));
        }
        // This existing bounded admission checks every region, not only trial
        // temperatures. Its immutable borrow ends before the model is published.
        model.prepare_caloric_branch(controls, &mut cancelled)?;
        let initial_h = initial.specific_enthalpy_j_kg();
        let evaluate = |t| -> Result<root::Sample, GasEnergyError> {
            let properties = model.thermodynamics()
                .evaluate_properties(Temperature::new(t), Pressure::new(pressure))
                .map_err(|e| failure("isobaric composition thermodynamics", e))?;
            Ok(root::Sample {
                residual: (properties.h().value() / properties.molar_mass().value() - initial_h) - q,
                tangent: properties.mass_specific_cp().value(),
                pressure_pa: pressure,
            })
        };
        let solved = root::solve(controls.minimum_temperature_k, controls.maximum_temperature_k,
            controls.energy_tolerance_j_kg, controls.maximum_iterations, evaluate, &mut cancelled)
            .map_err(map_root)?;
        poll(&mut cancelled)?;
        let point = updated_point(&old.point, solved.temperature_k, pressure)?;
        let final_state = model.evaluate(&point)?;
        let volume1 = mass / final_state.state().density;
        let work = pressure * (volume1 - volume0);
        let enthalpy = mass * (final_state.specific_enthalpy_j_kg() - initial_h);
        let stored = mass * (final_state.specific_internal_energy_j_kg()
            - initial.specific_internal_energy_j_kg());
        let enthalpy_residual = enthalpy - heat_j;
        let energy_residual = (stored + work) - heat_j;
        if !positive(volume1)
            || final_state.state().pressure.to_bits() != pressure.to_bits()
            || [work, enthalpy, stored, enthalpy_residual, energy_residual]
                .iter().any(|v| !v.is_finite())
            || (volume1 != volume0 && work == 0.0)
            || (heat_j != 0.0 && (enthalpy == 0.0
                || enthalpy.is_sign_positive() != heat_j.is_sign_positive()))
            || enthalpy_residual.abs() > tolerance_j || energy_residual.abs() > tolerance_j
        {
            return Err(GasEnergyError::Invalid(
                "isobaric composition endpoint did not close mass, enthalpy and actual volume work",
            ));
        }
        let report = GasIsobaricCompositionStep {
            pressure_pa: pressure,
            energy: GasCompositionStep {
                initial, final_state: final_state.clone(), mass_kg: mass,
                initial_volume_m3: volume0, volume_m3: volume1,
                supplied_heat_j: heat_j, work_by_gas_j: work,
                stored_energy_change_j: stored, energy_residual_j: energy_residual,
                iterations: solved.iterations,
            },
            enthalpy_change_j: enthalpy, enthalpy_residual_j: enthalpy_residual,
        };
        poll(&mut cancelled)?;
        self.parcel.model = model;
        self.parcel.point = point;
        self.parcel.state = final_state;
        self.parcel.volume_m3 = volume1;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
