//! Typed thermal material points for the original finite spatial owners.
//!
//! A thermal solve retains every constitutive history and freezes the actual
//! physical interval at each trial temperature. Inspection has no interval and
//! never re-integrates creep or recovery. The common thermal response preserves
//! physical heat separately from the variational flow potential and work defect.

use super::{
    FiniteJ2Response, FiniteJ2State, FiniteJ2ThermalError, FiniteJ2ThermalIncrement,
    FiniteJ2Thermodynamics, FiniteThermalJ2Material, FrozenThermalFiniteJ2,
};

/// One endpoint-temperature material point with complete, explicitly typed memory.
/// All deformation arguments are total F, and all energies use original volume.
pub trait FiniteThermalInelasticPoint {
    type State: Clone + core::fmt::Debug;

    /// Implicit response and analytic total-F tangent at frozen accepted memory.
    /// An inspection-only point must return DurationRequired for kinetic flow.
    fn response(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<FiniteJ2Response<Self::State>, FiniteJ2ThermalError>;

    /// Identical constitutive return without tangent assembly for line search.
    fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<([f64; 9], f64), FiniteJ2ThermalError>;

    /// Inspect the complete accepted state without flow, recovery or projection.
    fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2ThermalError>;

    /// Free energy, entropy and mechanical internal energy at fixed history.
    fn thermodynamics(
        &self,
        f: &[f64; 9],
        state: &Self::State,
    ) -> Result<FiniteJ2Thermodynamics, FiniteJ2ThermalError>;
}

/// A complete finite thermal family, retaining its actual history and kinetics.
///
/// Consumers freeze each physical step at its endpoint temperature and repeatedly
/// evaluate the same accepted history during Newton and thermal iteration. They
/// commit the entire returned State only when the surrounding step succeeds.
/// Caloric deposition is physical dissipation plus reversible exchange; the flow
/// potential and signed mechanical endpoint-work defect are separate accounts.
pub trait FiniteThermalInelasticMaterial {
    type State: Clone + core::fmt::Debug;
    type Frozen: FiniteThermalInelasticPoint<State = Self::State>;

    fn initial_state(&self) -> Self::State;
    fn reference_temperature_k(&self) -> f64;
    fn temperature_range_k(&self) -> [f64; 2];

    /// Add this conservative original-volume bound to the caloric capacity.
    /// A coupled owner must require their total to be strictly positive.
    fn minimum_mechanical_heat_capacity_j_m3_k(&self) -> f64;

    /// Inspect history and thermodynamics without inventing a physical interval.
    fn at_temperature_for_inspection(
        &self,
        temperature_k: f64,
    ) -> Result<Self::Frozen, FiniteJ2ThermalError>;

    /// Freeze the actual positive physical interval at the trial temperature.
    fn at_temperature_for_step(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<Self::Frozen, FiniteJ2ThermalError>;

    /// Check OLD history at OLD temperature, perform one native constitutive
    /// return, and retain the original thermal work/entropy/heat accounting.
    #[allow(clippy::too_many_arguments)]
    fn checked_increment_for_step(
        &self,
        old_f: &[f64; 9],
        old: &Self::State,
        old_temperature_k: f64,
        new_f: &[f64; 9],
        new_temperature_k: f64,
        yield_tolerance_pa: f64,
        duration_s: f64,
    ) -> Result<FiniteJ2ThermalIncrement<Self::State>, FiniteJ2ThermalError>;
}

// Existing native and source-admitted users keep their exact original owner.
impl<T: FiniteThermalJ2Material + ?Sized> FiniteThermalInelasticMaterial for T {
    type State = FiniteJ2State;
    type Frozen = FrozenThermalFiniteJ2;

    fn initial_state(&self) -> Self::State {
        self.thermal_finite_j2().initial_state()
    }

    fn reference_temperature_k(&self) -> f64 {
        self.thermal_finite_j2()
            .parameters()
            .reference_temperature_k
    }

    fn temperature_range_k(&self) -> [f64; 2] {
        self.thermal_finite_j2().parameters().temperature_range_k
    }

    fn minimum_mechanical_heat_capacity_j_m3_k(&self) -> f64 {
        self.thermal_finite_j2()
            .minimum_mechanical_heat_capacity_j_m3_k()
    }

    fn at_temperature_for_inspection(
        &self,
        temperature_k: f64,
    ) -> Result<Self::Frozen, FiniteJ2ThermalError> {
        self.thermal_finite_j2()
            .at_temperature_for_inspection(temperature_k)
    }

    fn at_temperature_for_step(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<Self::Frozen, FiniteJ2ThermalError> {
        self.thermal_finite_j2()
            .at_temperature_for_step(temperature_k, duration_s)
    }

    fn checked_increment_for_step(
        &self,
        old_f: &[f64; 9],
        old: &Self::State,
        old_temperature_k: f64,
        new_f: &[f64; 9],
        new_temperature_k: f64,
        yield_tolerance_pa: f64,
        duration_s: f64,
    ) -> Result<FiniteJ2ThermalIncrement<Self::State>, FiniteJ2ThermalError> {
        self.thermal_finite_j2().checked_increment_for_step(
            old_f,
            old,
            old_temperature_k,
            new_f,
            new_temperature_k,
            yield_tolerance_pa,
            duration_s,
        )
    }
}

impl FiniteThermalInelasticPoint for FrozenThermalFiniteJ2 {
    type State = FiniteJ2State;

    fn response(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<FiniteJ2Response<Self::State>, FiniteJ2ThermalError> {
        Self::response(self, f, old)
    }

    fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<([f64; 9], f64), FiniteJ2ThermalError> {
        Self::stress_and_potential(self, f, old)
    }

    fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2ThermalError> {
        Self::validate_committed(self, f, state, yield_tolerance_pa)
    }

    fn thermodynamics(
        &self,
        f: &[f64; 9],
        state: &Self::State,
    ) -> Result<FiniteJ2Thermodynamics, FiniteJ2ThermalError> {
        Self::thermodynamics(self, f, state)
    }
}
