//! Atomic gas-material replacement with one extensive energy target.
//!
//! This is the composition-change boundary for a uniform parcel. A caller may
//! supply a new thermochemical model with different mole fractions or a different
//! admitted species set. Mass remains fixed. The original replacement conserves
//! volume and internal energy; its heat/work variant reconstructs directly at
//! U_start + Q - W and the supplied final volume. Temperature and pressure come
//! from the NEW model at that target, then all transport and acoustic properties
//! resolve from the same state.
//!
//! This operation is deliberately not a mixing-rate, diffusion, valve, reaction,
//! or chemical-equilibrium model. The caller owns the physical process that chose
//! the new composition and must account for species mass flux/chemical work when
//! using it inside an open-system simulation.

use super::ClosedGasParcel;
use super::super::{
    GasEnergyError, GasEnergySolve, ResolvedThermochemicalGasState,
    ThermochemicalGasModel, poll, positive, updated_point,
};

/// Accepted replacement of the parcel's thermochemical/material model.
#[derive(Clone, Debug, PartialEq)]
pub struct GasCompositionStep {
    /// Complete source-resolved state before the replacement.
    pub initial: ResolvedThermochemicalGasState,
    /// Complete source-resolved state after the replacement.
    pub final_state: ResolvedThermochemicalGasState,
    /// Conserved parcel mass [kg].
    pub mass_kg: f64,
    /// Initial parcel volume [m3].
    pub initial_volume_m3: f64,
    /// Accepted parcel volume [m3]; unchanged in the original fixed-volume APIs.
    pub volume_m3: f64,
    /// Signed heat supplied while changing composition [J]; zero for the
    /// original conserved-energy replacement.
    pub supplied_heat_j: f64,
    /// Signed work delivered by the gas [J]; zero in the fixed-volume heat APIs.
    pub work_by_gas_j: f64,
    /// Actual total internal-energy change [J].
    pub stored_energy_change_j: f64,
    /// Actual stored change plus delivered work minus supplied heat [J].
    pub energy_residual_j: f64,
    /// Interior temperature-inversion iterations.
    pub iterations: usize,
}

impl ClosedGasParcel {
    /// Replace composition/species/source data while conserving mass, volume and
    /// total internal energy. The new model therefore determines the resulting
    /// temperature, pressure, viscosity, conductivity, sound speed and impedance.
    ///
    /// This is the data-only operation needed for scenarios such as replacing an
    /// air mixture with helium-enriched, CO2-enriched, or otherwise exotic gas
    /// material data without editing the mechanics/acoustics code. The supplied
    /// model may have a different admitted species set, provided its source laws
    /// support the conserved density/internal-energy state and retained external
    /// query axes.
    ///
    /// The old parcel is unchanged on source refusal, unbracketed energy,
    /// exhausted work, arithmetic failure or cancellation.
    pub fn replace_model_conserving_mass_energy(
        &mut self,
        model: ThermochemicalGasModel,
        controls: GasEnergySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasCompositionStep, GasEnergyError> {
        self.replace_model_with_heat(model, 0.0, controls, cancelled)
    }

    /// Replace the material model and admit a signed heat exchange in ONE
    /// original density/energy reconstruction. Mass and volume stay fixed.
    ///
    /// This avoids a fictitious intermediate at the new composition but old
    /// energy (or old composition but new energy). A coupled reaction/heat
    /// owner can therefore query its actual nonadiabatic endpoint, including
    /// formation references, without requiring an unsupported split state.
    /// The caller owns stoichiometry, kinetics and the equal/opposite thermal
    /// debit; this primitive supplies no additional chemical heat source.
    /// Refusal or cancellation preserves the complete accepted parcel.
    pub fn replace_model_with_heat(
        &mut self,
        model: ThermochemicalGasModel,
        heat_j: f64,
        controls: GasEnergySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasCompositionStep, GasEnergyError> {
        self.replace_model_with_heat_work_volume(
            model, self.volume_m3(), heat_j, 0.0, controls, cancelled,
        )
    }

    /// Reconstruct composition, volume, heat and work at ONE material endpoint.
    ///
    /// Mass stays fixed and the original caloric inversion targets U_start+Q-W.
    /// Positive W is work delivered by the gas. The enclosing mechanical owner
    /// must supply and independently close that work against its actual boundary
    /// motion; this material primitive does not infer a pressure quadrature,
    /// trajectory or entropy-production law. It creates no intermediate split
    /// composition/temperature/volume state and adds no separate formation heat.
    /// Failed source, energy, volume or cancellation admission changes nothing.
    #[allow(clippy::too_many_arguments)]
    pub fn replace_model_with_heat_work_volume(
        &mut self,
        model: ThermochemicalGasModel,
        volume: f64,
        heat_j: f64,
        work_by_gas_j: f64,
        controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasCompositionStep, GasEnergyError> {
        poll(&mut cancelled)?;
        let mass = self.parcel.mass_kg;
        let initial_volume = self.parcel.volume_m3;
        if !positive(mass) || !positive(volume) {
            return Err(GasEnergyError::Invalid("positive represented parcel mass and volume required"));
        }
        let initial = self.parcel.state.clone();
        let specific_u = initial.specific_internal_energy_j_kg();
        let exchange = heat_j - work_by_gas_j;
        let specific_exchange = exchange / mass;
        let target_u = specific_u + specific_exchange;
        if !specific_u.is_finite() || !heat_j.is_finite() || !work_by_gas_j.is_finite()
            || !exchange.is_finite() || !target_u.is_finite()
            || (heat_j != 0.0 && exchange == -work_by_gas_j)
            || (work_by_gas_j != 0.0 && exchange == heat_j)
            || (exchange != 0.0 && (specific_exchange == 0.0 || target_u == specific_u))
        {
            return Err(GasEnergyError::Invalid(
                "finite represented initial energy, heat and work increments required",
            ));
        }
        let density = mass / volume;
        if !positive(density)
            || (volume != initial_volume && density == mass / initial_volume)
        {
            return Err(GasEnergyError::Invalid("parcel endpoint density is unrepresentable"));
        }
        let solution = model.reconstruct_density_energy(
            &self.parcel.point,
            density,
            target_u,
            controls,
            &mut cancelled,
        )?;
        let final_state = solution.state().clone();
        let stored = mass * (final_state.specific_internal_energy_j_kg() - specific_u);
        let residual = (stored + work_by_gas_j) - heat_j;
        let tolerance_j = mass * controls.energy_tolerance_j_kg;
        if !stored.is_finite() || !residual.is_finite() || !positive(tolerance_j)
            || residual.abs() > tolerance_j
        {
            return Err(GasEnergyError::Invalid(
                "replacement model does not close represented parcel energy, heat and work",
            ));
        }
        let point = updated_point(
            &self.parcel.point,
            final_state.state().temperature,
            final_state.state().pressure,
        )?;
        let iterations = solution.iterations();
        poll(&mut cancelled)?;
        self.parcel.model = model;
        self.parcel.point = point;
        self.parcel.state = final_state.clone();
        self.parcel.volume_m3 = volume;
        Ok(GasCompositionStep {
            initial,
            final_state,
            mass_kg: mass,
            initial_volume_m3: initial_volume,
            volume_m3: volume,
            supplied_heat_j: heat_j,
            work_by_gas_j,
            stored_energy_change_j: stored,
            energy_residual_j: residual,
            iterations,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermochemical_gas::tests::{model, point};

    fn controls() -> GasEnergySolve {
        GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1500.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        }
    }

    fn parcel() -> ClosedGasParcel {
        ClosedGasParcel::new(
            model(vec![0.75, 0.25], 0.001),
            point(400.0, 100000.0).with("material-lot", 42.0).unwrap(),
            0.001,
        ).unwrap()
    }

    #[test]
    fn composition_swap_preserves_extensive_storage_and_changes_acoustics() {
        let mut gas = parcel();
        let before = gas.clone();
        let report = gas.replace_model_conserving_mass_energy(
            model(vec![0.25, 0.75], 0.001), controls(), || false,
        ).unwrap();
        assert_eq!(gas.mass_kg(), before.mass_kg());
        assert_eq!(gas.volume_m3(), before.volume_m3());
        assert_eq!(gas.point().axes().get("material-lot"), Some(&42.0));
        assert!(report.stored_energy_change_j.abs()
            <= gas.mass_kg() * controls().energy_tolerance_j_kg);
        assert_ne!(gas.state().state().pressure, before.state().state().pressure);
        assert_ne!(gas.state().state().sound_speed, before.state().state().sound_speed);
        assert_ne!(gas.state().state().dynamic_viscosity, before.state().state().dynamic_viscosity);
        assert_ne!(gas.state().state().thermal_conductivity, before.state().state().thermal_conductivity);
    }

    #[test]
    fn identical_model_is_storage_preserving_and_source_complete() {
        let mut gas = parcel();
        let before = gas.clone();
        let report = gas.replace_model_conserving_mass_energy(
            model(vec![0.75, 0.25], 0.001), controls(), || false,
        ).unwrap();
        // The new model re-solves T(rho, u), so the state agrees within the
        // inversion's declared energy tolerance, not bitwise: |dT| <= tol / cv
        // with the ideal-gas cv = R / (gamma - 1) of the accepted state.
        let (after, initial) = (gas.state().state(), before.state().state());
        let cv = initial.specific_gas_constant / (initial.gamma - 1.0);
        let temperature_tolerance = controls().energy_tolerance_j_kg / cv;
        assert!((after.temperature - initial.temperature).abs() <= temperature_tolerance);
        assert!(
            (after.pressure / initial.pressure - 1.0).abs()
                <= temperature_tolerance / initial.temperature + 72.0 * f64::EPSILON
        );
        assert!((after.density / initial.density - 1.0).abs() <= 64.0 * f64::EPSILON);
        assert_eq!(gas.mass_kg(), before.mass_kg());
        assert_eq!(gas.volume_m3(), before.volume_m3());
        assert!(report.stored_energy_change_j.abs()
            <= gas.mass_kg() * controls().energy_tolerance_j_kg);
        assert_eq!(gas.state().components().len(), 2);
    }

    #[test]
    fn failed_or_cancelled_replacement_is_atomic() {
        let mut gas = parcel();
        let before = gas.clone();
        // A bracket that provably excludes the replacement's actual solution.
        let solved = gas.clone().replace_model_conserving_mass_energy(
            model(vec![0.25, 0.75], 0.001), controls(), || false,
        ).unwrap().final_state.state().temperature;
        let narrow = GasEnergySolve { minimum_temperature_k: solved + 50.0,
            maximum_temperature_k: solved + 100.0, ..controls() };
        assert!(gas.replace_model_conserving_mass_energy(
            model(vec![0.25, 0.75], 0.001), narrow, || false,
        ).is_err());
        assert_eq!(gas, before);
        assert!(matches!(gas.replace_model_conserving_mass_energy(
            model(vec![0.25, 0.75], 0.001), controls(), || true,
        ), Err(GasEnergyError::Cancelled)));
        assert_eq!(gas, before);
    }

    #[test]
    fn finite_volume_heat_and_work_reconstruct_one_isothermal_material_endpoint() {
        let mut gas = parcel();
        let before = gas.clone();
        let ratio: f64 = 1.2;
        let initial = before.state().state();
        // Independent reversible ideal-gas isotherm: Q=W=m R T ln(V1/V0).
        // Total internal energy depends only on T for this frozen ideal mixture.
        let work = before.mass_kg() * initial.specific_gas_constant
            * initial.temperature * ratio.ln();
        let expected_entropy = work / initial.temperature;
        let report = gas.replace_model_with_heat_work_volume(
            before.model().clone(), before.volume_m3() * ratio,
            work, work, controls(), || false,
        ).unwrap();
        assert_eq!(gas.mass_kg(), before.mass_kg());
        assert_eq!(gas.volume_m3(), before.volume_m3() * ratio);
        assert_eq!(report.initial_volume_m3, before.volume_m3());
        assert_eq!(report.volume_m3, gas.volume_m3());
        assert_eq!(report.supplied_heat_j, work);
        assert_eq!(report.work_by_gas_j, work);
        assert!(report.stored_energy_change_j.abs()
            <= before.mass_kg() * controls().energy_tolerance_j_kg);
        assert!(report.energy_residual_j.abs()
            <= before.mass_kg() * controls().energy_tolerance_j_kg);
        assert!((gas.state().state().temperature - initial.temperature).abs() < 1e-8);
        assert!((gas.state().state().pressure / (initial.pressure / ratio) - 1.0).abs()
            < 1e-10);
        assert!((gas.state().state().density / (initial.density / ratio) - 1.0).abs()
            < 1e-12);
        let entropy = gas.mass_kg() * (
            gas.state().specific_entropy_j_kg_k().unwrap()
            - before.state().specific_entropy_j_kg_k().unwrap()
        );
        assert!((entropy - expected_entropy).abs() < 1e-9);
        assert_eq!(gas.point().axes().get("material-lot"), Some(&42.0));
    }

    #[test]
    fn failed_moving_material_endpoint_preserves_volume_sources_and_energy() {
        let mut gas = parcel();
        let before = gas.clone();
        assert!(gas.replace_model_with_heat_work_volume(
            before.model().clone(), before.volume_m3() * 1.2,
            1e12, 0.0, controls(), || false,
        ).is_err());
        assert_eq!(gas, before);
        let mut polls = 0;
        let cancelled = gas.replace_model_with_heat_work_volume(
            before.model().clone(), before.volume_m3() * 1.2,
            0.0, 0.0, controls(), || { polls += 1; polls >= 4 },
        );
        assert!(matches!(cancelled, Err(GasEnergyError::Cancelled)));
        assert!(polls >= 4);
        assert_eq!(gas, before);
    }
}
