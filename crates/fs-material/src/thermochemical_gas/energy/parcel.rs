//! Closed gas storage with material-driven heating and reversible boundary work.
//!
//! One existing isochoric parcel owns mass, geometry and the resolved state.
//! Heating holds its current volume; a separate isentropic operation changes
//! volume without heat. Composition, source selection and extra query axes
//! survive both operations. The separate [`transfer`] operation exchanges mass
//! between volumes; no piston inertia or combustion is implied. The opt-in [`midpoint`] update couples heat and volume work in
//! conservative coordinates; it is not an exact finite isentropic operation.
//! The [`isobaric`] path heats against a quasistatic constant-pressure load,
//! with enthalpy and the actual pressure-volume work checked independently.

/// Conservative composition/model replacement at fixed mass, volume and total internal energy.
pub mod composition;
/// Constant-pressure heating with independent enthalpy and boundary-work balances.
pub mod isobaric;
/// Work-conjugate mechanical loads using the common acoustic compressibility.
pub mod mechanics;
/// Simultaneous heating and volume work through the original material laws.
pub mod midpoint;
/// Conservative paired gas mass, species and enthalpy transfer at fixed volumes.
pub mod transfer;

use super::entropy::{GasEntropyError, GasEntropySolve, map_context};
use super::{GasEnergyError, GasEnergySolve, GasHeatStep, IsochoricGasParcel,
    ResolvedThermochemicalGasState, ThermochemicalGasModel, poll, positive, updated_point};
use fs_matdb::QueryPoint;

/// Closed, variable-volume gas. Storage and heating retain the existing owner;
/// no second species registry, caloric law, or acoustic medium is introduced.
#[derive(Clone, Debug, PartialEq)]
pub struct ClosedGasParcel {
    parcel: IsochoricGasParcel,
}

/// Accepted reversible adiabatic volume change. Positive work leaves the gas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasVolumeStep {
    /// Initial volume [m3].
    pub initial_volume_m3: f64,
    /// Accepted final volume [m3].
    pub final_volume_m3: f64,
    /// Original-law `m Delta u` [J], canceling common formation references
    /// before rounding rather than subtracting absolute stored energies.
    pub stored_energy_change_j: f64,
    /// Work delivered by the gas [J], equal to minus its stored-energy change.
    /// This is the thermodynamic work on the admitted isentrope, NOT an
    /// independent pressure-quadrature or measured mechanical-work certificate.
    pub work_by_gas_j: f64,
    /// Actual `m*(s_final-s_initial)` [J/K], bounded by the caller's tolerance.
    pub entropy_change_j_k: f64,
    /// Interior entropy-inversion iterations, excluding bounded endpoint checks.
    pub iterations: usize,
}
impl ClosedGasParcel {
    /// Resolve initial T/p and obtain mass from that same material density.
    pub fn new(model: ThermochemicalGasModel, point: QueryPoint, volume_m3: f64)
        -> Result<Self, GasEnergyError>
    {
        IsochoricGasParcel::new(model, point, volume_m3).map(Self::from_isochoric)
    }
    /// Preserve an already resolved parcel without resampling or losing sources.
    #[must_use]
    pub fn from_isochoric(parcel: IsochoricGasParcel) -> Self { Self { parcel } }
    /// Freeze the current geometry without re-evaluating any material property.
    #[must_use]
    pub fn into_isochoric(self) -> IsochoricGasParcel { self.parcel }
    /// Common gas state for momentum, thermal and acoustic consumers.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalGasState { self.parcel.state() }
    /// Conserved mass [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 { self.parcel.mass_kg() }
    /// Current volume [m3].
    #[must_use]
    pub const fn volume_m3(&self) -> f64 { self.parcel.volume_m3() }
    /// Updated T/p and all retained external material conditions.
    #[must_use]
    pub const fn point(&self) -> &QueryPoint { self.parcel.point() }
    /// Original active species, caloric/transport sources and composition.
    /// Immutable access lets reaction/flow consumers validate common laws
    /// without introducing a second gas state or changing accepted storage.
    #[must_use]
    pub const fn model(&self) -> &ThermochemicalGasModel { &self.parcel.model }
    /// Prescribed signed heat at the CURRENT held volume. The existing energy
    /// closure, source-domain checks and atomic publication are unchanged.
    pub fn add_heat(&mut self, heat_j: f64, controls: GasEnergySolve, cancelled: impl FnMut() -> bool)
        -> Result<GasHeatStep, GasEnergyError>
    {
        self.parcel.add_heat(heat_j, controls, cancelled)
    }

    /// Reconstruct a frozen-composition endpoint from volume and actual Delta U.
    ///
    /// This is a conservative-coordinate adapter to the ORIGINAL temperature-
    /// anchored caloric inversion, not a heat/work or mechanical integrator.
    /// The enclosing physical owner must close the returned actual energy change
    /// against its heat and boundary work. Common formation references cancel
    /// before evaluating either the inversion or the returned joule increment.
    ///
    /// Exact accepted mass, species, source model and auxiliary query axes are
    /// retained. The original parcel is immutable, including on cancellation or
    /// any unsupported density/temperature/source/energy trial.
    pub fn at_volume_internal_energy_increment(
        &self, volume_m3: f64, energy_change_j: f64, controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, f64), GasEnergyError> {
        poll(&mut cancelled)?;
        let old = &self.parcel;
        let density = old.mass_kg / volume_m3;
        let requested = energy_change_j / old.mass_kg;
        let tolerance = old.mass_kg * controls.energy_tolerance_j_kg;
        if !positive(volume_m3) || !positive(density) || !positive(tolerance)
            || !energy_change_j.is_finite() || !requested.is_finite()
            || (energy_change_j != 0.0 && requested == 0.0)
            || (volume_m3 != old.volume_m3 && density == old.mass_kg / old.volume_m3)
        {
            return Err(GasEnergyError::Invalid(
                "positive represented volume/density and finite energy coordinate required",
            ));
        }
        let solved = old.model.reconstruct_density_internal_energy_increment(
            &old.point, density, requested, controls, &mut cancelled,
        )?;
        let state = solved.state();
        let stored = old.mass_kg * old.model.thermodynamics()
            .specific_internal_energy_change(
                fs_qty::Temperature::new(old.state.state().temperature),
                fs_qty::Temperature::new(state.state().temperature),
            ).map_err(|error| crate::thermochemical_gas::failure(
                "frozen volume/energy coordinate", error,
            ))?.value();
        if !stored.is_finite() || (stored - energy_change_j).abs() > tolerance {
            return Err(GasEnergyError::Invalid(
                "frozen endpoint does not retain the original caloric increment",
            ));
        }
        let mut parcel = old.clone();
        parcel.point = updated_point(&old.point, state.state().temperature, state.state().pressure)?;
        parcel.state = state.clone();
        parcel.volume_m3 = volume_m3;
        poll(&mut cancelled)?;
        Ok((Self::from_isochoric(parcel), stored))
    }

    /// Reversible adiabatic compression/expansion at conserved mass and entropy.
    /// The caller must admit one continuous, stable gas branch over the supplied
    /// temperature bracket, INCLUDING the initial state. The endpoint is fully
    /// re-resolved through the same material law; no constant-gamma rescaling
    /// or extrapolated source property is substituted.
    ///
    /// Heat and volume changes are separate physical operations. Sequentially
    /// composing them does not establish a time integrator for simultaneous
    /// heating and motion. Work here excludes any external-pressure load, piston
    /// inertia or friction; the mechanical consumer must account for those.
    ///
    /// Failed source queries, unrepresentable work, exhausted iteration budgets
    /// and cancellation preserve the entire parcel, including its geometry.
    pub fn set_volume_isentropic(
        &mut self, volume_m3: f64, controls: GasEntropySolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasVolumeStep, GasEntropyError> {
        controls.validate()?;
        poll(&mut cancelled).map_err(map_context)?;
        let old = &self.parcel;
        let temperature = old.state.state().temperature;
        if !positive(volume_m3) || temperature < controls.minimum_temperature_k
            || temperature > controls.maximum_temperature_k
        {
            return Err(GasEntropyError::Invalid("positive finite volume and a bracket containing the initial temperature required"));
        }
        let tolerance = old.mass_kg * controls.entropy_tolerance_j_kg_k;
        if !positive(tolerance) {
            return Err(GasEntropyError::Invalid("total entropy tolerance is unrepresentable"));
        }
        let mut report = GasVolumeStep {
            initial_volume_m3: old.volume_m3, final_volume_m3: volume_m3,
            stored_energy_change_j: 0.0, work_by_gas_j: 0.0,
            entropy_change_j_k: 0.0, iterations: 0,
        };
        if volume_m3 == old.volume_m3 {
            poll(&mut cancelled).map_err(map_context)?;
            return Ok(report);
        }
        let density = old.mass_kg / volume_m3;
        if !positive(density) || density == old.mass_kg / old.volume_m3 {
            return Err(GasEntropyError::Invalid("volume change does not give a representable positive density change"));
        }
        let solved = old.model.reconstruct_density_entropy_change(&old.point, density,
            0.0, controls, &mut cancelled)?;
        let stored = old.mass_kg * old.model.thermodynamics().specific_internal_energy_change(
            fs_qty::Temperature::new(temperature),
            fs_qty::Temperature::new(solved.state().state().temperature),
        ).map_err(|error| crate::thermochemical_gas::failure("isentropic volume work", error))?.value();
        let entropy_change = old.mass_kg * solved.residual_j_kg_k();
        // Positive cv and pressure imply heating on compression and cooling on
        // expansion. Do not accept vanished work under an overly loose tolerance.
        if !stored.is_finite() || !entropy_change.is_finite() || entropy_change.abs() > tolerance
            || (volume_m3 < old.volume_m3 && stored <= 0.0)
            || (volume_m3 > old.volume_m3 && stored >= 0.0)
        {
            return Err(GasEntropyError::Invalid("isentropic volume work or entropy balance is unrepresentable/inconsistent"));
        }
        let point = updated_point(&old.point, solved.state().state().temperature,
            solved.state().state().pressure).map_err(map_context)?;
        report.stored_energy_change_j = stored;
        report.work_by_gas_j = -stored;
        report.entropy_change_j_k = entropy_change;
        report.iterations = solved.iterations();
        poll(&mut cancelled).map_err(map_context)?;
        self.parcel.state = solved.into_state();
        self.parcel.point = point;
        self.parcel.volume_m3 = volume_m3;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;


#[cfg(test)]
mod conservative_coordinate_tests {
    use super::*;
    use crate::thermochemical_gas::tests::{model, point};

    #[test]
    fn frozen_volume_coordinate_retains_exact_mass_and_isothermal_source_state() {
        let gas = ClosedGasParcel::new(
            model(vec![0.75, 0.25], 0.001),
            point(400.0, 100000.0).with("lot", 42.0).unwrap(), 0.001,
        ).unwrap();
        let solve = GasEnergySolve { minimum_temperature_k: 200.0,
            maximum_temperature_k: 1500.0, energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128 };
        // Independent ideal-gas isotherm: zero Delta U, p proportional to 1/V.
        let (candidate, change) = gas.at_volume_internal_energy_increment(
            0.0012, 0.0, solve, || false,
        ).unwrap();
        assert_eq!(candidate.mass_kg().to_bits(), gas.mass_kg().to_bits());
        assert_eq!(candidate.model(), gas.model());
        assert_eq!(candidate.point().axes().get("lot"), Some(&42.0));
        assert_eq!(candidate.volume_m3(), 0.0012);
        assert!(change.abs() <= gas.mass_kg() * solve.energy_tolerance_j_kg);
        assert!((candidate.state().state().temperature - 400.0).abs() < 1e-8);
        assert!((candidate.state().state().pressure - 100000.0 / 1.2).abs() < 1e-6);
        let before = gas.clone();
        assert!(gas.at_volume_internal_energy_increment(-1.0, 0.0, solve, || false).is_err());
        assert!(gas.at_volume_internal_energy_increment(0.0012, 1e20, solve, || false).is_err());
        assert!(matches!(gas.at_volume_internal_energy_increment(
            0.0012, 0.0, solve, || true,
        ), Err(GasEnergyError::Cancelled)));
        assert_eq!(gas, before);
    }
}
