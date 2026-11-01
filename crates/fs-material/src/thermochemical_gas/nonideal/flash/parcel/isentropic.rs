//! Reversible adiabatic volume changes of a closed phase-changing mixture.
//!
//! The original entropy/volume owner resolves pressure, temperature and phase
//! amounts at each prescribed volume. Mass and the complete feed are unchanged.
//! Boundary work BY the mixture is inferred as minus the original source Delta U
//! under the reversible, adiabatic, quasistatic assumptions. It is not a second,
//! independently integrated work balance. There is no piston inertia, heat input,
//! phase kinetics, slip, reaction, frozen gamma or invented mixture heat capacity.

use super::{
    E, ThermochemicalPengRobinsonModel, admit_inventory, poll, positive, source_specific_volume,
};
use crate::thermochemical_gas::nonideal::flash::enthalpy::spacing;
use crate::thermochemical_gas::nonideal::flash::equilibrium::{
    ResolvedThermochemicalPhaseState, entropy::volume::PhaseEntropyVolumeControls,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::eos::peng_robinson::phase_screen::PengRobinsonPhaseScreenOptionsV1;

/// One accepted quasistatic volume change with original source caloric changes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsentropicMixtureVolumeStep {
    /// Prescribed final physical volume minus the accepted initial volume [m3].
    pub volume_change_m3: f64,
    /// Original source-integrated internal-energy change [J].
    pub stored_energy_change_j: f64,
    /// Independently source-integrated enthalpy change [J].
    pub enthalpy_change_j: f64,
    /// Reversible adiabatic work BY the mixture, defined as minus Delta U [J].
    /// Positive expansion exports energy. This is not independent work quadrature.
    pub boundary_work_j: f64,
    /// Original source Delta S, also the zero-entropy-constraint residual [J/K].
    /// Numerical residuals remain visible and are never clipped to zero.
    pub entropy_change_j_k: f64,
    /// Original EOS total volume minus the prescribed published volume [m3].
    pub volume_residual_m3: f64,
    /// Actual attempted pressure inversions in the original S,V solve.
    pub pressure_evaluations: usize,
    /// Actual original PT selections, including the reference and all inner work.
    pub pt_evaluations: usize,
    /// Original screen work plus split ceilings; not an observed EOS-call count.
    pub eos_evaluation_upper_bound: usize,
}

/// Conserved original mixture inventory following prescribed isentropic volumes.
///
/// The accepted physical volume is the caller's target, never silently replaced
/// by the nearby EOS volume. Source residuals remain explicit. Phase compositions
/// and amounts may change while mass and total species composition stay fixed.
/// Cloning retains the full immutable source and accepted state for exact retry.
#[derive(Clone, Debug, PartialEq)]
pub struct IsentropicMixtureEquilibriumParcel {
    model: ThermochemicalPengRobinsonModel,
    state: ResolvedThermochemicalPhaseState,
    mass_kg: f64,
    volume_m3: f64,
    specific_volume_m3_kg: f64,
}

impl IsentropicMixtureEquilibriumParcel {
    /// Resolve the initial original PT state and its matching physical volume.
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
            return Err(E::Invalid(
                "positive finite isentropic mixture mass required",
            ));
        }
        let state = model.evaluate_phase_equilibrium(
            temperature,
            pressure,
            phase_options,
            &mut cancelled,
        )?;
        let volume_m3 = mass_kg * source_specific_volume(&state);
        let parcel = Self::checked(model, state, mass_kg, volume_m3)?;
        poll(&mut cancelled)?;
        Ok(parcel)
    }

    /// Resolve initial mass, volume and original-reference specific entropy.
    /// When a large absolute entropy datum cannot resolve the tolerance, use PT
    /// initialization followed by original source-increment volume changes.
    pub fn from_entropy(
        model: ThermochemicalPengRobinsonModel,
        mass_kg: f64,
        volume_m3: f64,
        specific_entropy_j_kg_k: f64,
        controls: PhaseEntropyVolumeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, E> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let specific_volume = admit_inventory(mass_kg, volume_m3)?;
        let tolerance_m3 = mass_kg * controls.specific_volume_tolerance_m3_kg;
        if !positive(tolerance_m3) || spacing(volume_m3) > tolerance_m3 {
            return Err(E::Invalid(
                "isentropic mixture volume tolerance is unrepresentable",
            ));
        }
        let solved = model.phase_equilibrium_entropy_volume(
            specific_volume,
            specific_entropy_j_kg_k,
            controls,
            &mut cancelled,
        )?;
        let residual = mass_kg * source_specific_volume(solved.state()) - volume_m3;
        if !residual.is_finite() || residual.abs() > tolerance_m3 {
            return Err(E::Invalid(
                "isentropic mixture initial volume balance refused",
            ));
        }
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
        admit_stored_state(&state, mass_kg)?;
        Ok(Self {
            model,
            state,
            mass_kg,
            volume_m3,
            specific_volume_m3_kg,
        })
    }

    /// Immutable original species, reference conventions and interaction laws.
    #[must_use]
    pub const fn model(&self) -> &ThermochemicalPengRobinsonModel {
        &self.model
    }
    /// Complete accepted phase state, original source receipts and finite screen.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalPhaseState {
        &self.state
    }
    /// Conserved total feed mass [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Accepted prescribed physical volume [m3], within the reported EOS residual.
    #[must_use]
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Physical inventory density [kg/m3], with source density retained separately.
    #[must_use]
    pub fn density_kg_m3(&self) -> f64 {
        self.mass_kg / self.volume_m3
    }
    /// Total internal energy on its original reference [J]. For small changes
    /// under large datums, use the original-source increment in the step report.
    #[must_use]
    pub fn internal_energy_j(&self) -> f64 {
        self.mass_kg * self.state.specific_internal_energy_j_kg()
    }
    /// Total enthalpy on its original reference [J].
    #[must_use]
    pub fn enthalpy_j(&self) -> f64 {
        self.mass_kg * self.state.specific_enthalpy_j_kg()
    }
    /// Total entropy on its original species reference [J/K].
    #[must_use]
    pub fn entropy_j_k(&self) -> f64 {
        self.mass_kg * self.state.specific_entropy_j_kg_k()
    }

    /// Prescribe a new volume through reversible adiabatic equilibrium.
    ///
    /// State and physical volume publish together only after source, entropy,
    /// volume, representability and work-direction checks. Small nonzero motion
    /// reserves both entropy and volume resolution against its pressure-work
    /// scale. These safeguards are not certified relative work-error bounds.
    /// Motion below the retained initial volume resolution refuses.
    ///
    /// Exact unchanged volume preserves every accepted bit and performs no phase
    /// solve; it reports the retained volume residual without refining that state.
    /// Nonzero motion requires the accepted original phase policy, so reevaluating
    /// its reference cannot silently introduce an unreported caloric change.
    /// An exact rest accepts any valid controls without applying a new policy.
    /// No duration, piston force, acceleration, heat-transfer law or kinetics is
    /// inferred from this sequence of quasistatic thermodynamic endpoints.
    pub fn change_volume(
        &mut self,
        final_volume_m3: f64,
        mut controls: PhaseEntropyVolumeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsentropicMixtureVolumeStep, E> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let final_specific_volume = admit_inventory(self.mass_kg, final_volume_m3)?;
        let volume_change_m3 = final_volume_m3 - self.volume_m3;
        let initial_volume_residual =
            self.mass_kg * source_specific_volume(&self.state) - self.volume_m3;
        if volume_change_m3 == 0.0 {
            let report = IsentropicMixtureVolumeStep {
                volume_change_m3: 0.0,
                stored_energy_change_j: 0.0,
                enthalpy_change_j: 0.0,
                boundary_work_j: 0.0,
                entropy_change_j_k: 0.0,
                volume_residual_m3: initial_volume_residual,
                pressure_evaluations: 0,
                pt_evaluations: 0,
                eos_evaluation_upper_bound: 0,
            };
            poll(&mut cancelled)?;
            return Ok(report);
        }
        if controls.phase.phase_options() != self.state.phase().options() {
            return Err(E::Invalid(
                "continuing isentropic mixture requires its accepted original phase policy",
            ));
        }
        let change_m3_kg = final_specific_volume - self.specific_volume_m3_kg;
        let min_pressure = self
            .state
            .pressure_pa()
            .min(controls.pressure_bracket[0].value());
        let max_pressure = self
            .state
            .pressure_pa()
            .max(controls.pressure_bracket[1].value());
        let max_temperature = self
            .state
            .temperature_k()
            .max(controls.phase.temperature_bracket[1].value());
        let work_scale_j_kg = min_pressure * change_m3_kg.abs();
        let entropy_allowance = 0.125 * work_scale_j_kg / max_temperature;
        let volume_allowance = 0.125 * work_scale_j_kg / max_pressure;
        if change_m3_kg == 0.0
            || !positive(work_scale_j_kg)
            || !positive(entropy_allowance)
            || !positive(volume_allowance)
            || initial_volume_residual.abs() > 0.125 * volume_change_m3.abs()
        {
            return Err(E::Invalid(
                "isentropic mixture motion is below retained source resolution",
            ));
        }
        controls.phase.entropy_tolerance_j_kg_k = controls
            .phase
            .entropy_tolerance_j_kg_k
            .min(entropy_allowance);
        controls.specific_volume_tolerance_m3_kg = controls
            .specific_volume_tolerance_m3_kg
            .min(0.25 * change_m3_kg.abs())
            .min(volume_allowance);
        let tolerance_j_k = self.mass_kg * controls.phase.entropy_tolerance_j_kg_k;
        let tolerance_m3 = self.mass_kg * controls.specific_volume_tolerance_m3_kg;
        if !positive(tolerance_j_k)
            || !positive(tolerance_m3)
            || spacing(final_volume_m3) > tolerance_m3
        {
            return Err(E::Invalid(
                "isentropic mixture motion tolerances are unrepresentable",
            ));
        }
        let solved = self.model.phase_equilibrium_entropy_volume_increment(
            Temperature::new(self.state.temperature_k()),
            Pressure::new(self.state.pressure_pa()),
            final_specific_volume,
            0.0,
            controls,
            &mut cancelled,
        )?;
        let delta_u = solved.internal_energy_change_j_kg().ok_or(E::Invalid(
            "isentropic mixture source internal-energy increment is missing",
        ))?;
        let delta_h = solved.enthalpy_change_j_kg().ok_or(E::Invalid(
            "isentropic mixture source enthalpy increment is missing",
        ))?;
        let delta_s = solved.accepted_increment_j_kg_k().ok_or(E::Invalid(
            "isentropic mixture source entropy increment is missing",
        ))?;
        let stored_energy_change_j = self.mass_kg * delta_u;
        let enthalpy_change_j = self.mass_kg * delta_h;
        let boundary_work_j = -stored_energy_change_j;
        let entropy_change_j_k = self.mass_kg * delta_s;
        let volume_residual_m3 =
            self.mass_kg * source_specific_volume(solved.state()) - final_volume_m3;
        if !stored_energy_change_j.is_finite()
            || !enthalpy_change_j.is_finite()
            || !entropy_change_j_k.is_finite()
            || !volume_residual_m3.is_finite()
            || boundary_work_j == 0.0
            || boundary_work_j.signum() != volume_change_m3.signum()
            || entropy_change_j_k.abs() > tolerance_j_k
            || volume_residual_m3.abs() > tolerance_m3
            || solved.state().constitutive_law() != self.state.constitutive_law()
        {
            return Err(E::Invalid(
                "isentropic mixture entropy, volume or work direction refused",
            ));
        }
        admit_stored_state(solved.state(), self.mass_kg)?;
        let report = IsentropicMixtureVolumeStep {
            volume_change_m3,
            stored_energy_change_j,
            enthalpy_change_j,
            boundary_work_j,
            entropy_change_j_k,
            volume_residual_m3,
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

fn admit_stored_state(state: &ResolvedThermochemicalPhaseState, mass_kg: f64) -> Result<(), E> {
    if !positive(mass_kg * source_specific_volume(state))
        || [
            state.specific_internal_energy_j_kg(),
            state.specific_enthalpy_j_kg(),
            state.specific_entropy_j_kg_k(),
        ]
        .into_iter()
        .any(|value| !(mass_kg * value).is_finite())
    {
        return Err(E::Invalid(
            "isentropic mixture extensive state is unrepresentable",
        ));
    }
    Ok(())
}
