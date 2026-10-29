//! Finite source-backed surrounding gas on the existing live string step.
//!
//! Material viscous loss and reversible thermoelastic exchange heat the string;
//! aerodynamic drag heats the gas.
//! The existing NASA-9 parcel owns caloric inversion, EOS and transport. The
//! acoustic owner then rebuilds drag and its compact observer from that same
//! accepted gas, preserving modal coordinates and the previous acceleration.
//! This is a first-order partitioned, well-mixed fixed-volume model. The
//! `exchange` module adds reciprocal source-resolved thermal transfer. Neither
//! path introduces a chamber boundary, radiation loading or gas pressure work.
//! The `flow` module composes the original conservative species/enthalpy nozzle
//! owner with that continuing acoustic state when gas enters or leaves.

/// Ideal electrical source work on the same finite-gas/string transaction.
pub mod electrical;
/// Reciprocal source-resolved conduction on the same finite-gas/string step.
pub mod exchange;
/// Source-resolved pressure-driven species exchange on the continuing string.
pub mod flow;
/// Source-resolved endpoint chemistry on the continuing string and finite gas.
pub mod reactive;

use super::{
    AcousticRealizeError, ThermalMaterialStringFrame, ThermalMaterialStringRuntime,
    ThermalStringError, linear_accelerations,
};
use fs_material::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, GasHeatStep, IsochoricGasParcel,
};

/// Two independent prescribed thermal sources during one audio sample.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FiniteGasStringHeat {
    /// Signed heat into the solid [J], excluding its computed viscous loss.
    pub string_heat_j: f64,
    /// Signed heat into the gas [J], excluding computed aerodynamic drag.
    pub gas_heat_j: f64,
}

/// Source inversion controls and a physical partition-step limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasStringSolve {
    /// Existing source-temperature bracket, absolute energy tolerance and cap.
    pub gas: GasEnergySolve,
    /// Maximum absolute gas temperature change in one accepted sample [K].
    /// Limits the frozen-initial-medium approximation, not an accuracy proof.
    pub maximum_gas_temperature_increment_k: f64,
}

/// No error publishes either the string or its finite surrounding gas.
#[derive(Debug)]
pub enum FiniteGasStringError {
    /// Incompatible state, source accounting or coupling limit.
    Admission(&'static str),
    /// Original string/material/acoustic refusal.
    String(ThermalStringError),
    /// Original gas/source inversion refusal.
    Gas(GasEnergyError),
    /// Original chemical mechanism, source, nonlinear or consumption refusal.
    Reaction(crate::cavity_phs::piston::reactive::GasReactionError),
}
impl core::fmt::Display for FiniteGasStringError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-gas material string: {self:?}")
    }
}
impl std::error::Error for FiniteGasStringError {}
impl From<ThermalStringError> for FiniteGasStringError {
    fn from(error: ThermalStringError) -> Self {
        Self::String(error)
    }
}
impl From<GasEnergyError> for FiniteGasStringError {
    fn from(error: GasEnergyError) -> Self {
        Self::Gas(error)
    }
}
impl From<crate::cavity_phs::piston::reactive::GasReactionError> for FiniteGasStringError {
    fn from(error: crate::cavity_phs::piston::reactive::GasReactionError) -> Self {
        Self::Reaction(error)
    }
}

/// Actual mechanical, solid-thermal and gas storage at one accepted epoch.
#[derive(Debug)]
pub struct FiniteGasStringFrame {
    /// Existing solid/string account. Its air loss is the INTERNAL transfer
    /// `drag_heat_j` in the combined account below, not a second energy debit.
    pub string: ThermalMaterialStringFrame,
    /// Existing gas-owner account, including drag, prescribed gas heat and any
    /// separately reported reciprocal solid-gas transfer.
    pub gas: GasHeatStep,
    /// Aerodynamic loss deposited exactly once into the gas [J].
    pub drag_heat_j: f64,
    /// External gas heat only [J], excluding that internal transfer.
    pub external_gas_heat_j: f64,
    /// Work + two prescribed heats - complete solid/mechanical/gas storage changes - authored loss
    /// - the explicitly reported damping-partition roundoff residual [J].
    /// All storage changes are measured from accepted endpoint states.
    pub energy_balance_residual_j: f64,
    /// Sum of original absolute solve allowances and arithmetic roundoff [J].
    pub energy_tolerance_j: f64,
}

struct FiniteGasStringCandidate {
    runtime: ThermalMaterialStringRuntime,
    gas: IsochoricGasParcel,
    string: ThermalMaterialStringFrame,
    gas_step: GasHeatStep,
}

fn check(cx: &fs_exec::Cx<'_>) -> Result<(), FiniteGasStringError> {
    cx.checkpoint().map_err(|_| {
        FiniteGasStringError::String(ThermalStringError::Acoustic(
            AcousticRealizeError::Cancelled,
        ))
    })
}

fn modal_energy(runtime: &ThermalMaterialStringRuntime) -> f64 {
    runtime
        .mechanical
        .modes()
        .iter()
        .zip(runtime.mechanical.states())
        .map(|(mode, state)| {
            0.5 * (state.velocity_m_sqrt_kg_per_s.powi(2)
                + (mode.angular_frequency_rad_s * state.displacement_m_sqrt_kg).powi(2))
        })
        .sum()
}

impl ThermalMaterialStringRuntime {
    /// Advance vibration, solid heating and finite gas heating atomically.
    ///
    /// The parcel's complete initial snapshot must equal the runtime's ambient.
    /// Initialize with `LinearMaterialStringRuntime::try_new_with_gas`, or
    /// explicitly `rebind_gas` before starting this coupled evolution.
    /// Species, source receipts, mass and volume remain owned by the parcel.
    ///
    /// One exact mechanical sample freezes the initial material viscosity and
    /// gas drag. Material loss enters the solid once; air drag enters the gas
    /// once. Their existing energy owners reconstruct both temperatures and
    /// the final gas EOS/transport. Final drag, pressure observation and retained
    /// acceleration use that same endpoint gas at ONE new epoch. Rejected
    /// gas/source, temperature, observer, energy or cancellation checks preserve
    /// both accepted states and consume neither source input nor audio sample.
    ///
    /// The finite gas is an explicitly well-mixed local thermal inventory;
    /// its volume is not a resonator geometry. This model omits conductive/
    /// convective exchange between string and gas, gas motion/pressure work,
    /// radiation loading and propagation of a changing medium. The compact
    /// pressure observer remains diagnostic and is not deposited as extra heat.
    /// The first-order split requires timestep refinement; temperature bounds
    /// alone do not establish thermodynamic or acoustic accuracy.
    pub fn step_with_isochoric_gas(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        controls: FiniteGasStringSolve,
    ) -> Result<FiniteGasStringFrame, FiniteGasStringError> {
        self.admit_finite_gas_step(cx, expected_epoch, gas, heat, controls)?;
        let mut candidate = self.clone();
        let mut next_gas = gas.clone();
        let string = candidate.step(cx, expected_epoch, forces, heat.string_heat_j)?;
        let drag_heat_j = string.vibration.dissipation.air_loss_j;
        let total_gas_heat = heat.gas_heat_j + drag_heat_j;
        if !total_gas_heat.is_finite()
            || (drag_heat_j != 0.0 && total_gas_heat == heat.gas_heat_j)
            || (heat.gas_heat_j != 0.0 && total_gas_heat == drag_heat_j)
        {
            return Err(FiniteGasStringError::Admission(
                "gas source and aerodynamic heat overflow or disappear when combined",
            ));
        }
        let gas_step =
            next_gas.add_heat(total_gas_heat, controls.gas, || cx.checkpoint().is_err())?;
        self.publish_finite_gas(
            cx,
            forces,
            gas,
            heat,
            controls,
            FiniteGasStringCandidate {
                runtime: candidate,
                gas: next_gas,
                string,
                gas_step,
            },
        )
    }

    fn admit_finite_gas_step(
        &self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        gas: &IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        controls: FiniteGasStringSolve,
    ) -> Result<(), FiniteGasStringError> {
        check(cx)?;
        if expected_epoch != self.mechanical.epoch() {
            return Err(ThermalStringError::Epoch {
                expected: expected_epoch,
                actual: self.mechanical.epoch(),
            }
            .into());
        }
        controls.gas.validate()?;
        if !heat.string_heat_j.is_finite()
            || !heat.gas_heat_j.is_finite()
            || !controls.maximum_gas_temperature_increment_k.is_finite()
            || controls.maximum_gas_temperature_increment_k <= 0.0
            || self.mechanical.ambient() != gas.state().state()
        {
            return Err(FiniteGasStringError::Admission(
                "finite heats, a positive temperature-increment limit and one identical initial gas snapshot required",
            ));
        }
        Ok(())
    }

    fn publish_finite_gas(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        controls: FiniteGasStringSolve,
        proposal: FiniteGasStringCandidate,
    ) -> Result<FiniteGasStringFrame, FiniteGasStringError> {
        let FiniteGasStringCandidate {
            runtime: mut candidate,
            gas: next_gas,
            mut string,
            gas_step,
        } = proposal;
        let initial_internal_energy = gas.state().specific_internal_energy_j_kg();
        let initial_mechanical_energy = modal_energy(self);
        let drag_heat_j = string.vibration.dissipation.air_loss_j;
        if (next_gas.state().state().temperature - gas.state().state().temperature).abs()
            > controls.maximum_gas_temperature_increment_k
        {
            return Err(FiniteGasStringError::Admission(
                "finite gas exceeded its per-sample temperature-increment limit",
            ));
        }

        self.refresh_finite_gas_medium(
            cx,
            forces,
            &mut candidate,
            &mut string,
            next_gas.state().state(),
        )?;

        let final_internal_energy = next_gas.state().specific_internal_energy_j_kg();
        let gas_storage = gas.mass_kg() * (final_internal_energy - initial_internal_energy);
        let solid_storage = self.mechanical.specimen().mass_kg()
            * (candidate.thermal.specific_enthalpy_j_kg() - self.thermal.specific_enthalpy_j_kg());
        let mechanical_storage = modal_energy(&candidate) - initial_mechanical_energy
            + string.axial_free_energy_change_j
            + string.mechanical_entropy_energy_change_j;
        let dissipation = &string.vibration.dissipation;
        let input = string.vibration.acoustic.input_work_j + heat.string_heat_j + heat.gas_heat_j;
        let residual = input
            - (mechanical_storage
                + solid_storage
                + gas_storage
                + dissipation.authored_loss_j
                + dissipation.roundoff_residual_j);
        let gas_tolerance = gas.mass_kg() * controls.gas.energy_tolerance_j_kg;
        let scale = input.abs()
            + heat.gas_heat_j.abs()
            + drag_heat_j
            + mechanical_storage.abs()
            + solid_storage.abs()
            + gas_storage.abs()
            + gas.mass_kg() * (initial_internal_energy.abs() + final_internal_energy.abs());
        let tolerance = string.energy_roundoff_tolerance_j
            + string.thermal_solve_tolerance_j
            + gas_tolerance
            + 32.0 * f64::EPSILON * scale;
        if !residual.is_finite()
            || !tolerance.is_finite()
            || gas_tolerance <= 0.0
            || residual.abs() > tolerance
        {
            return Err(FiniteGasStringError::Admission(
                "independently measured string/solid/gas energy balance did not close",
            ));
        }
        check(cx)?;
        // No fallible operation separates publication. Destructors run after
        // both exclusive owners have been replaced.
        let _old_string = core::mem::replace(self, candidate);
        let _old_gas = core::mem::replace(gas, next_gas);
        Ok(FiniteGasStringFrame {
            string,
            gas: gas_step,
            drag_heat_j,
            external_gas_heat_j: heat.gas_heat_j,
            energy_balance_residual_j: residual,
            energy_tolerance_j: tolerance,
        })
    }

    // The surrounding transaction may heat AND exchange species before it has
    // its final medium. Always reconstruct from the original sample history,
    // never from an intermediate observer evaluation in the same sample.
    fn refresh_finite_gas_medium(
        &self,
        cx: &fs_exec::Cx<'_>,
        forces: &[f64],
        candidate: &mut Self,
        string: &mut ThermalMaterialStringFrame,
        gas: &fs_material::gas::GasState,
    ) -> Result<(), FiniteGasStringError> {
        // Rebuild at this sample's publication epoch; no extra physical instant.
        // The rebind retains q and p and owes zero elastic parameter work.
        candidate.mechanical.epoch = self.mechanical.epoch();
        candidate
            .mechanical
            .rebind_gas(cx, gas)
            .map_err(ThermalStringError::Acoustic)?;
        candidate.mechanical.member.observer.previous_acceleration = self
            .mechanical
            .member
            .observer
            .previous_acceleration
            .clone();
        let initial = candidate
            .mechanical
            .member
            .observer
            .previous_acceleration
            .is_none()
            .then(|| linear_accelerations(&self.mechanical.member.model, forces));
        string.vibration.acoustic.observer_pressure_pa = candidate
            .mechanical
            .member
            .observer
            .observe(
                linear_accelerations(&candidate.mechanical.member.model, forces),
                initial,
                candidate.mechanical.member.model.sample_period_s(),
            )
            .map_err(ThermalStringError::Acoustic)?;
        Ok(())
    }
}
