//! Adiabatic well-mixed withdrawal from a rigid phase-equilibrium inventory.
//!
//! Every infinitesimal outgoing portion has the current bulk composition and
//! specific entropy. With no heat or shaft work, the remaining specific entropy
//! stays constant and its specific volume becomes V/m. The original S,V owner
//! supplies that endpoint, including any redistribution between phases.
//! Integrating d(m*u) = h*dm gives the total exported enthalpy from the original
//! internal-energy loss; freezing the initial h over a finite portion is avoided.
//!
//! This assumes bulk entrainment of both equilibrium phases, no slip, and
//! negligible stored kinetic/potential energy. It supplies no duration, port
//! rate, valve/nozzle loss, phase-selective venting or downstream equilibrium.

use super::{
    E, IsochoricMixtureEquilibriumParcel, admit_inventory, poll, positive, source_specific_volume,
};
use crate::thermochemical_gas::nonideal::flash::enthalpy::spacing;
use crate::thermochemical_gas::nonideal::flash::equilibrium::entropy::volume::PhaseEntropyVolumeControls;
use fs_qty::{Pressure, Temperature};

/// Original entropy/volume policy and the physical mass-roundoff allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureEquilibriumDischargeControl {
    /// Original source, phase, entropy/volume tolerances and nested work caps.
    pub equilibrium: PhaseEntropyVolumeControls,
    /// Maximum absolute removed-minus-requested mass residual [kg].
    pub mass_tolerance_kg: f64,
}
impl MixtureEquilibriumDischargeControl {
    /// Use an explicit original S,V bracket with a 1e-12 kg mass allowance.
    #[must_use]
    pub fn new(equilibrium: PhaseEntropyVolumeControls) -> Self {
        Self {
            equilibrium,
            mass_tolerance_kg: 1e-12,
        }
    }

    /// Validate the complete original solve and mass allowances.
    pub fn validate(self) -> Result<(), E> {
        self.equilibrium.validate()?;
        if !positive(self.mass_tolerance_kg) {
            return Err(E::Invalid(
                "positive finite discharge mass tolerance required",
            ));
        }
        Ok(())
    }
}

/// One accepted fixed-volume withdrawal and its original caloric accounting.
///
/// Absolute exported energy/entropy retain the species reference convention.
/// The specific source changes and exported-enthalpy adjustment cancel those
/// datums before arithmetic. Mean enthalpy is an integral over the withdrawal,
/// not the enthalpy of one asserted outlet temperature/pressure/phase state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureEquilibriumDischargeStep {
    /// Caller-prescribed outgoing mass [kg].
    pub requested_mass_kg: f64,
    /// Original mass minus the representable accepted remaining mass [kg].
    pub removed_mass_kg: f64,
    /// Removed mass minus requested mass [kg], never silently clipped.
    pub mass_residual_kg: f64,
    /// Original source-integrated specific internal-energy change of storage [J/kg].
    pub remaining_specific_internal_energy_change_j_kg: f64,
    /// Original source-integrated specific enthalpy change of storage [J/kg].
    pub remaining_specific_enthalpy_change_j_kg: f64,
    /// Total outgoing enthalpy minus removed_mass * original inlet h [J].
    /// This source-datum-independent correction includes cooling during withdrawal.
    pub exported_enthalpy_adjustment_j: f64,
    /// Integrated total outgoing enthalpy [J], on the original species reference.
    /// It need not be positive under an arbitrary formation-energy convention.
    pub exported_enthalpy_j: f64,
    /// Negative of exported enthalpy [J], inferred from the adiabatic first law.
    /// This is not an independent quadrature or residual.
    pub stored_energy_change_j: f64,
    /// Portion-mean outgoing h [J/kg]; absent for exact zero withdrawal.
    pub mean_exported_enthalpy_j_kg: Option<f64>,
    /// Portion-mean h minus the original inlet h [J/kg], datum independent.
    pub mean_exported_enthalpy_change_j_kg: Option<f64>,
    /// Removed mass times the original specific entropy [J/K].
    /// Entropy generation outside the vessel is not assigned to this receipt.
    pub exported_entropy_j_k: f64,
    /// Remaining mass times the original source Delta-s [J/K].
    /// The ideal well-mixed path is isentropic; its numerical residual stays visible.
    pub entropy_balance_residual_j_k: f64,
    /// Original final EOS volume minus the unchanged physical vessel volume [m³].
    pub volume_residual_m3: f64,
    /// Actual original pressure/entropy inversions.
    pub pressure_evaluations: usize,
    /// Actual original PT selections, including the reference and all inner solves.
    pub pt_evaluations: usize,
    /// Original EOS work ceiling, not an observed number of cubic evaluations.
    pub eos_evaluation_upper_bound: usize,
}

impl IsochoricMixtureEquilibriumParcel {
    /// Withdraw a prescribed mass through an adiabatic well-mixed equilibrium path.
    ///
    /// The complete original feed composition is retained by both the vessel and
    /// every outgoing portion. Final mass, phase state and density publish together
    /// after the native S,V closure, mass/entropy/volume checks and final cancellation.
    /// Physical vessel volume never changes and no boundary work is added.
    ///
    /// Positive downstream pressure must not exceed either accepted endpoint.
    /// These endpoint checks are necessary, not a certification of the continuous
    /// pressure path or a resolved pressure-driven flow rate.
    ///
    /// Rounding of the prescribed mass is explicit and bounded; a portion that
    /// cannot change the retained inventory/volume coordinates is refused. Exact
    /// zero withdrawal preserves all state bits and performs no source evaluation.
    pub fn discharge_adiabatic(
        &mut self,
        requested_mass_kg: f64,
        downstream_pressure: Pressure,
        controls: MixtureEquilibriumDischargeControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MixtureEquilibriumDischargeStep, E> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !requested_mass_kg.is_finite()
            || requested_mass_kg < 0.0
            || !positive(downstream_pressure.value())
        {
            return Err(E::Invalid(
                "finite nonnegative discharge and positive downstream pressure required",
            ));
        }
        let initial_source_volume = source_specific_volume(&self.state);
        let initial_volume_residual = self.mass_kg * initial_source_volume - self.volume_m3;
        if requested_mass_kg == 0.0 {
            let report = MixtureEquilibriumDischargeStep {
                requested_mass_kg,
                removed_mass_kg: 0.0,
                mass_residual_kg: 0.0,
                remaining_specific_internal_energy_change_j_kg: 0.0,
                remaining_specific_enthalpy_change_j_kg: 0.0,
                exported_enthalpy_adjustment_j: 0.0,
                exported_enthalpy_j: 0.0,
                stored_energy_change_j: 0.0,
                mean_exported_enthalpy_j_kg: None,
                mean_exported_enthalpy_change_j_kg: None,
                exported_entropy_j_k: 0.0,
                entropy_balance_residual_j_k: 0.0,
                volume_residual_m3: initial_volume_residual,
                pressure_evaluations: 0,
                pt_evaluations: 0,
                eos_evaluation_upper_bound: 0,
            };
            poll(&mut cancelled)?;
            return Ok(report);
        }
        let mut equilibrium = controls.equilibrium;
        if equilibrium.phase.phase_options() != self.state.phase().options() {
            return Err(E::Invalid(
                "discharge requires the accepted original phase policy",
            ));
        }
        if downstream_pressure.value() >= self.state.pressure_pa() {
            return Err(E::Invalid(
                "discharge requires positive initial pressure head",
            ));
        }
        let remaining_mass = self.mass_kg - requested_mass_kg;
        let removed_mass = self.mass_kg - remaining_mass;
        let mass_residual = removed_mass - requested_mass_kg;
        if !positive(remaining_mass)
            || remaining_mass >= self.mass_kg
            || !positive(removed_mass)
            || removed_mass >= self.mass_kg
            || mass_residual.abs() > controls.mass_tolerance_kg.min(0.125 * requested_mass_kg)
            || spacing(self.mass_kg) > controls.mass_tolerance_kg
        {
            return Err(E::Invalid(
                "discharge must leave a resolved positive inventory",
            ));
        }
        let final_specific_volume = admit_inventory(remaining_mass, self.volume_m3)?;
        let volume_change = final_specific_volume - self.specific_volume_m3_kg;
        let min_pressure = self
            .state
            .pressure_pa()
            .min(equilibrium.pressure_bracket[0].value());
        let max_pressure = self
            .state
            .pressure_pa()
            .max(equilibrium.pressure_bracket[1].value());
        let max_temperature = self
            .state
            .temperature_k()
            .max(equilibrium.phase.temperature_bracket[1].value());
        let work_scale = min_pressure * volume_change;
        let entropy_allowance = 0.125 * work_scale / max_temperature;
        let volume_allowance = 0.125 * work_scale / max_pressure;
        // Reserve resolution against the finite portion's p*Delta-v scale.
        // These safeguards are not a certified relative exported-energy bound.
        if !positive(volume_change)
            || !positive(work_scale)
            || !positive(entropy_allowance)
            || !positive(volume_allowance)
            || (initial_source_volume - self.specific_volume_m3_kg).abs() > 0.125 * volume_change
        {
            return Err(E::Invalid(
                "discharge is below retained original-state resolution",
            ));
        }
        equilibrium.phase.entropy_tolerance_j_kg_k = equilibrium
            .phase
            .entropy_tolerance_j_kg_k
            .min(entropy_allowance);
        equilibrium.specific_volume_tolerance_m3_kg = equilibrium
            .specific_volume_tolerance_m3_kg
            .min(0.25 * volume_change)
            .min(volume_allowance);
        let entropy_tolerance = remaining_mass * equilibrium.phase.entropy_tolerance_j_kg_k;
        let volume_tolerance = remaining_mass * equilibrium.specific_volume_tolerance_m3_kg;
        if !positive(entropy_tolerance)
            || !positive(volume_tolerance)
            || spacing(self.volume_m3) > volume_tolerance
        {
            return Err(E::Invalid(
                "discharge conservation tolerances are unrepresentable",
            ));
        }
        let solved = self.model.phase_equilibrium_entropy_volume_increment(
            Temperature::new(self.state.temperature_k()),
            Pressure::new(self.state.pressure_pa()),
            final_specific_volume,
            0.0,
            equilibrium,
            &mut cancelled,
        )?;
        let delta_u = solved
            .internal_energy_change_j_kg()
            .ok_or(E::Invalid("discharge source Delta-u is missing"))?;
        let delta_h = solved
            .enthalpy_change_j_kg()
            .ok_or(E::Invalid("discharge source Delta-h is missing"))?;
        let delta_s = solved
            .accepted_increment_j_kg_k()
            .ok_or(E::Invalid("discharge source Delta-s is missing"))?;
        let entropy_residual = remaining_mass * delta_s;
        let volume_residual =
            remaining_mass * source_specific_volume(solved.state()) - self.volume_m3;
        // Integral H_out = dm*u0 - m1*Delta-u. Retain the departure from dm*h0
        // separately so formation datums never enter this physically important
        // finite-withdrawal correction.
        let initial_flow_work = self.state.pressure_pa() * initial_source_volume;
        let adjustment = (-remaining_mass).mul_add(delta_u, -removed_mass * initial_flow_work);
        let mean_change = adjustment / removed_mass;
        let mean_enthalpy = self.state.specific_enthalpy_j_kg() + mean_change;
        let exported_enthalpy =
            removed_mass.mul_add(self.state.specific_enthalpy_j_kg(), adjustment);
        let exported_entropy = removed_mass * self.state.specific_entropy_j_kg_k();
        if [
            delta_u,
            delta_h,
            delta_s,
            entropy_residual,
            volume_residual,
            adjustment,
            mean_change,
            mean_enthalpy,
            exported_enthalpy,
            exported_entropy,
            remaining_mass * solved.state().specific_internal_energy_j_kg(),
        ]
        .into_iter()
        .any(|value| !value.is_finite())
            || delta_u >= 0.0
            || entropy_residual.abs() > entropy_tolerance
            || volume_residual.abs() > volume_tolerance
            || solved.state().pressure_pa() > self.state.pressure_pa()
            || solved.state().pressure_pa() < downstream_pressure.value()
            || solved.state().constitutive_law() != self.state.constitutive_law()
        {
            return Err(E::Invalid(
                "discharge source, entropy, volume or pressure endpoint refused",
            ));
        }
        let report = MixtureEquilibriumDischargeStep {
            requested_mass_kg,
            removed_mass_kg: removed_mass,
            mass_residual_kg: mass_residual,
            remaining_specific_internal_energy_change_j_kg: delta_u,
            remaining_specific_enthalpy_change_j_kg: delta_h,
            exported_enthalpy_adjustment_j: adjustment,
            exported_enthalpy_j: exported_enthalpy,
            stored_energy_change_j: -exported_enthalpy,
            mean_exported_enthalpy_j_kg: Some(mean_enthalpy),
            mean_exported_enthalpy_change_j_kg: Some(mean_change),
            exported_entropy_j_k: exported_entropy,
            entropy_balance_residual_j_k: entropy_residual,
            volume_residual_m3: volume_residual,
            pressure_evaluations: solved.pressure_evaluations(),
            pt_evaluations: solved.pt_evaluations(),
            eos_evaluation_upper_bound: solved.eos_evaluation_upper_bound(),
        };
        poll(&mut cancelled)?;
        self.state = solved.into_state();
        self.mass_kg = remaining_mass;
        self.specific_volume_m3_kg = final_specific_volume;
        Ok(report)
    }
}
