//! First-law-coupled volume changes of a closed phase-changing mixture.
//!
//! The caller supplies the final physical volume, heat INTO the mixture and
//! boundary work BY the mixture. Their difference is the required original
//! source internal-energy increment, `Delta U = Q - W_by`. The existing U,V
//! equilibrium owner then recovers pressure, temperature and phase amounts.
//! This module deliberately does not invent a pressure path or approximate
//! `integral p dV`: a piston, structural or flow owner must provide that work.
//! Mass and the complete species inventory remain fixed.

use super::{
    E, ThermochemicalPengRobinsonModel, admit_inventory, poll, positive, source_specific_volume,
};
use crate::thermochemical_gas::nonideal::flash::enthalpy::spacing;
use crate::thermochemical_gas::nonideal::flash::equilibrium::{
    ResolvedThermochemicalPhaseState, volume::PhaseEnergyVolumeControls,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::eos::peng_robinson::phase_screen::PengRobinsonPhaseScreenOptionsV1;

/// One accepted coupled heat/work/volume endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DrivenMixtureEquilibriumStep {
    /// Prescribed final physical volume minus the previously accepted volume [m3].
    pub volume_change_m3: f64,
    /// Heat supplied TO the mixture [J]. Positive heats the inventory.
    pub supplied_heat_j: f64,
    /// Mechanical boundary work done BY the mixture [J]. Positive expansion work
    /// removes energy from the inventory.
    pub boundary_work_by_mixture_j: f64,
    /// Requested first-law internal-energy change, Q - W_by [J].
    pub requested_internal_energy_change_j: f64,
    /// Original source-integrated internal-energy change [J].
    pub stored_energy_change_j: f64,
    /// `Delta U - (Q - W_by)` [J].
    pub first_law_residual_j: f64,
    /// Original EOS volume minus the prescribed physical volume [m3].
    pub volume_residual_m3: f64,
    /// Original source entropy change between equilibrium endpoints [J/K].
    ///
    /// This is a state difference, not an entropy-production audit: reservoir
    /// temperature and mechanical dissipation are intentionally outside this API.
    pub entropy_change_j_k: f64,
    /// Actual outer pressure inversions.
    pub pressure_evaluations: usize,
    /// Actual original PT phase selections, including the reference.
    pub pt_evaluations: usize,
    /// Conservative original EOS-work upper bound.
    pub eos_evaluation_upper_bound: usize,
}

/// Closed mixture inventory driven by externally accounted heat and boundary work.
///
/// This is the general conserved-energy counterpart of the specialized isochoric
/// heat parcel and isentropic volume parcel. The accepted physical volume is the
/// caller's value; the nearby EOS volume remains visible as a residual. No
/// pressure trajectory, piston dynamics, heat-transfer law, reaction or phase
/// kinetics is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct DrivenMixtureEquilibriumParcel {
    model: ThermochemicalPengRobinsonModel,
    state: ResolvedThermochemicalPhaseState,
    mass_kg: f64,
    volume_m3: f64,
    specific_volume_m3_kg: f64,
}

impl DrivenMixtureEquilibriumParcel {
    /// Resolve an initial original PT equilibrium and adopt its EOS volume.
    pub fn from_pressure_temperature(
        model: ThermochemicalPengRobinsonModel,
        mass_kg: f64,
        temperature: Temperature,
        pressure: Pressure,
        phase_options: PengRobinsonPhaseScreenOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, E> {
        poll(&mut cancelled)?;
        if !positive(mass_kg) {
            return Err(E::Invalid("positive finite driven-mixture mass required"));
        }
        let state =
            model.evaluate_phase_equilibrium(temperature, pressure, phase_options, &mut cancelled)?;
        let volume_m3 = mass_kg * source_specific_volume(&state);
        let parcel = Self::checked(model, state, mass_kg, volume_m3)?;
        poll(&mut cancelled)?;
        Ok(parcel)
    }

    /// Recover an initial state from absolute U,V on the original source datum.
    pub fn from_internal_energy(
        model: ThermochemicalPengRobinsonModel,
        mass_kg: f64,
        volume_m3: f64,
        specific_internal_energy_j_kg: f64,
        controls: PhaseEnergyVolumeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, E> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let specific_volume = admit_inventory(mass_kg, volume_m3)?;
        let solved = model.phase_equilibrium_internal_energy_volume(
            specific_volume,
            specific_internal_energy_j_kg,
            controls,
            &mut cancelled,
        )?;
        let parcel = Self::checked(model, solved.into_state(), mass_kg, volume_m3)?;
        poll(&mut cancelled)?;
        Ok(parcel)
    }

    fn checked(
        model: ThermochemicalPengRobinsonModel,
        state: ResolvedThermochemicalPhaseState,
        mass_kg: f64,
        volume_m3: f64,
    ) -> Result<Self, E> {
        let specific_volume_m3_kg = admit_inventory(mass_kg, volume_m3)?;
        if !(mass_kg * state.specific_internal_energy_j_kg()).is_finite()
            || !positive(mass_kg * source_specific_volume(&state))
        {
            return Err(E::Invalid("driven-mixture extensive state is unrepresentable"));
        }
        Ok(Self {
            model,
            state,
            mass_kg,
            volume_m3,
            specific_volume_m3_kg,
        })
    }

    #[must_use]
    pub const fn model(&self) -> &ThermochemicalPengRobinsonModel {
        &self.model
    }
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalPhaseState {
        &self.state
    }
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    #[must_use]
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    #[must_use]
    pub fn density_kg_m3(&self) -> f64 {
        self.mass_kg / self.volume_m3
    }
    #[must_use]
    pub fn internal_energy_j(&self) -> f64 {
        self.mass_kg * self.state.specific_internal_energy_j_kg()
    }

    /// Apply one externally accounted heat/work/volume change atomically.
    ///
    /// `heat_j` is positive INTO the mixture and `boundary_work_by_mixture_j`
    /// is positive OUT of it. The material owner therefore receives exactly
    /// `(heat - work)/mass` as its U increment. The complete candidate is
    /// published only after independent extensive first-law and volume checks.
    ///
    /// Exact no-op inputs preserve every accepted bit and perform no phase solve.
    /// Nonzero changes retain the accepted original phase policy. Failure, work
    /// exhaustion and final cancellation leave both state and physical volume
    /// unchanged.
    pub fn advance(
        &mut self,
        final_volume_m3: f64,
        heat_j: f64,
        boundary_work_by_mixture_j: f64,
        controls: PhaseEnergyVolumeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenMixtureEquilibriumStep, E> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !heat_j.is_finite() || !boundary_work_by_mixture_j.is_finite() {
            return Err(E::Invalid("finite driven-mixture heat and boundary work required"));
        }
        let final_specific_volume = admit_inventory(self.mass_kg, final_volume_m3)?;
        let volume_change_m3 = final_volume_m3 - self.volume_m3;
        let requested_internal_energy_change_j = heat_j - boundary_work_by_mixture_j;
        if !requested_internal_energy_change_j.is_finite() {
            return Err(E::Invalid("driven-mixture first-law target overflow"));
        }
        if volume_change_m3 == 0.0 && heat_j == 0.0 && boundary_work_by_mixture_j == 0.0 {
            let report = DrivenMixtureEquilibriumStep {
                volume_change_m3: 0.0,
                supplied_heat_j: 0.0,
                boundary_work_by_mixture_j: 0.0,
                requested_internal_energy_change_j: 0.0,
                stored_energy_change_j: 0.0,
                first_law_residual_j: 0.0,
                volume_residual_m3: self.mass_kg
                    * (source_specific_volume(&self.state) - self.specific_volume_m3_kg),
                entropy_change_j_k: 0.0,
                pressure_evaluations: 0,
                pt_evaluations: 0,
                eos_evaluation_upper_bound: 0,
            };
            poll(&mut cancelled)?;
            return Ok(report);
        }
        if controls.phase.phase_options() != self.state.phase().options() {
            return Err(E::Invalid(
                "continuing driven mixture requires its accepted original phase policy",
            ));
        }

        let change_j_kg = requested_internal_energy_change_j / self.mass_kg;
        let tolerance_j = self.mass_kg * controls.internal_energy_tolerance_j_kg;
        let tolerance_m3 = self.mass_kg * controls.specific_volume_tolerance_m3_kg;
        if !change_j_kg.is_finite()
            || !positive(tolerance_j)
            || !positive(tolerance_m3)
            || spacing(final_volume_m3) > tolerance_m3
            || requested_internal_energy_change_j != 0.0 && change_j_kg == 0.0
        {
            return Err(E::Invalid("driven-mixture target or tolerances are unrepresentable"));
        }

        let solved = self.model.phase_equilibrium_internal_energy_volume_increment(
            Temperature::new(self.state.temperature_k()),
            Pressure::new(self.state.pressure_pa()),
            final_specific_volume,
            change_j_kg,
            controls,
            &mut cancelled,
        )?;
        let increment = solved.accepted_increment_j_kg().ok_or(E::Invalid(
            "driven-mixture source internal-energy increment is missing",
        ))?;
        let entropy = solved.entropy_change_j_kg_k().ok_or(E::Invalid(
            "driven-mixture source entropy increment is missing",
        ))?;

        let stored_energy_change_j = self.mass_kg * increment;
        let first_law_residual_j =
            stored_energy_change_j - requested_internal_energy_change_j;
        let volume_residual_m3 =
            self.mass_kg * source_specific_volume(solved.state()) - final_volume_m3;
        let entropy_change_j_k = self.mass_kg * entropy;
        if !stored_energy_change_j.is_finite()
            || !first_law_residual_j.is_finite()
            || !volume_residual_m3.is_finite()
            || !entropy_change_j_k.is_finite()
            || first_law_residual_j.abs() > tolerance_j
            || volume_residual_m3.abs() > tolerance_m3
            || solved.state().constitutive_law() != self.state.constitutive_law()
            || !(self.mass_kg * solved.state().specific_internal_energy_j_kg()).is_finite()
        {
            return Err(E::Invalid(
                "driven-mixture first-law or volume balance refused",
            ));
        }

        let report = DrivenMixtureEquilibriumStep {
            volume_change_m3,
            supplied_heat_j: heat_j,
            boundary_work_by_mixture_j,
            requested_internal_energy_change_j,
            stored_energy_change_j,
            first_law_residual_j,
            volume_residual_m3,
            entropy_change_j_k,
            pressure_evaluations: solved.pressure_evaluations(),
            pt_evaluations: solved.pt_evaluations(),
            eos_evaluation_upper_bound: solved.eos_evaluation_upper_bound(),
        };
        poll(&mut cancelled)?;
        self.state = solved.into_state();
        self.volume_m3 = final_volume_m3;
        self.specific_volume_m3_kg = final_specific_volume;
        Ok(report)
    }
}
