//! Heat-driven mixture equilibrium in a closed, rigid vessel.
//!
//! The original source-bound U,V owner determines pressure, temperature and
//! phase split together. Mass, composition and vessel volume are fixed. Heat
//! changes the original total internal energy; latent energy is not an extra
//! source. This is instantaneous phase equilibrium, not evaporation kinetics,
//! wall transport, mechanical work, slip or an acoustic closure.

use super::equilibrium::ResolvedThermochemicalPhaseState;
use super::equilibrium::volume::PhaseEnergyVolumeControls;
use super::{E, ThermochemicalPengRobinsonModel, poll, positive};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::eos::peng_robinson::phase_screen::PengRobinsonPhaseScreenOptionsV1;

/// Closed phase-changing inventories against a constant pressure load.
pub mod isobaric;

/// One accepted closed-vessel heat increment and its independent source balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureEquilibriumHeatStep {
    /// Signed heat supplied to the vessel [J].
    pub supplied_heat_j: f64,
    /// Original source-integrated internal-energy increment [J].
    pub stored_energy_change_j: f64,
    /// Stored increment minus supplied heat [J].
    pub energy_residual_j: f64,
    /// Source-predicted volume minus the unchanged vessel volume [m3].
    pub volume_residual_m3: f64,
    /// Source-integrated entropy change [J/K], with species datums cancelled.
    /// Heat alone does not specify the reservoir temperature or entropy produced.
    pub entropy_change_j_k: f64,
    /// Original bounded pressure samples in this U,V solve.
    pub pressure_evaluations: usize,
    /// Original bounded PT phase selections, including the initial state.
    pub pt_evaluations: usize,
    /// Original screen work plus split work ceilings; not an observed EOS count.
    pub eos_evaluation_upper_bound: usize,
}

/// A fixed mixture inventory and volume with one accepted thermodynamic state.
///
/// The retained law includes every species amount and interaction parameter.
/// Heating may change phase compositions and phase amounts, but never the total
/// feed. The source-predicted volume is admitted to the explicit U,V tolerance;
/// the vessel volume is never changed to hide a solver residual.
#[derive(Clone, Debug, PartialEq)]
pub struct IsochoricMixtureEquilibriumParcel {
    model: ThermochemicalPengRobinsonModel,
    state: ResolvedThermochemicalPhaseState,
    mass_kg: f64,
    volume_m3: f64,
    specific_volume_m3_kg: f64,
}

impl IsochoricMixtureEquilibriumParcel {
    /// Resolve the initial PT phase state and define its fixed vessel volume.
    /// This avoids requiring a caller to subtract large formation-energy datums
    /// merely to initialize a physical temperature/pressure state.
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
            return Err(E::Invalid("positive finite mixture parcel mass required"));
        }
        let state = model.evaluate_phase_equilibrium(
            temperature,
            pressure,
            phase_options,
            &mut cancelled,
        )?;
        let specific_volume_m3_kg = source_specific_volume(&state);
        let volume_m3 = mass_kg * specific_volume_m3_kg;
        let parcel = Self::checked(model, state, mass_kg, volume_m3)?;
        poll(&mut cancelled)?;
        Ok(parcel)
    }

    /// Resolve the initial state from conserved mass, volume and specific U.
    /// Absolute U uses the original species reference convention. If that datum
    /// cannot resolve the requested tolerance, initialize with PT and then use
    /// source-integrated heat increments instead.
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
            return Err(E::Invalid("mixture parcel total energy is unrepresentable"));
        }
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
    /// Complete accepted phase state, source receipts and finite stability screen.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalPhaseState {
        &self.state
    }
    /// Conserved total feed mass [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Fixed physical vessel volume [m3], never adjusted by a solve.
    #[must_use]
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Conserved vessel density [kg/m3]. The source state differs only within
    /// the explicitly admitted volume residual; it is not overwritten.
    #[must_use]
    pub fn density_kg_m3(&self) -> f64 {
        self.mass_kg / self.volume_m3
    }
    /// Accepted total internal energy on the original reference [J]. For small
    /// changes under large datums, use the step's integrated increment receipt.
    #[must_use]
    pub fn internal_energy_j(&self) -> f64 {
        self.mass_kg * self.state.specific_internal_energy_j_kg()
    }

    /// Add signed heat and atomically publish the complete new phase state.
    ///
    /// Stable source increments cancel species formation and entropy datums
    /// before rounding. The U,V owner closes both constraints; this consumer
    /// independently checks its physical joule/volume accounts before mutation.
    /// Failure, exhausted budgets or late cancellation leave the parcel intact.
    /// A zero input performs no phase solve and preserves every stored bit.
    pub fn add_heat(
        &mut self,
        heat_j: f64,
        controls: PhaseEnergyVolumeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MixtureEquilibriumHeatStep, E> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !heat_j.is_finite() {
            return Err(E::Invalid("finite signed mixture parcel heat required"));
        }
        if heat_j == 0.0 {
            let report = MixtureEquilibriumHeatStep {
                supplied_heat_j: heat_j,
                stored_energy_change_j: 0.0,
                energy_residual_j: 0.0,
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
        let change_j_kg = heat_j / self.mass_kg;
        let tolerance_j = self.mass_kg * controls.internal_energy_tolerance_j_kg;
        let tolerance_m3 = self.mass_kg * controls.specific_volume_tolerance_m3_kg;
        if !change_j_kg.is_finite()
            || change_j_kg == 0.0
            || !positive(tolerance_j)
            || !positive(tolerance_m3)
        {
            return Err(E::Invalid(
                "mixture parcel heat or tolerances are unrepresentable",
            ));
        }
        let solved = self
            .model
            .phase_equilibrium_internal_energy_volume_increment(
                Temperature::new(self.state.temperature_k()),
                Pressure::new(self.state.pressure_pa()),
                self.specific_volume_m3_kg,
                change_j_kg,
                controls,
                &mut cancelled,
            )?;
        let increment = solved.accepted_increment_j_kg().ok_or(E::Invalid(
            "mixture parcel source energy increment receipt is missing",
        ))?;
        let entropy = solved.entropy_change_j_kg_k().ok_or(E::Invalid(
            "mixture parcel source entropy increment receipt is missing",
        ))?;
        let stored_energy_change_j = self.mass_kg * increment;
        let energy_residual_j = stored_energy_change_j - heat_j;
        let volume_residual_m3 =
            self.mass_kg * (source_specific_volume(solved.state()) - self.specific_volume_m3_kg);
        let entropy_change_j_k = self.mass_kg * entropy;
        if !stored_energy_change_j.is_finite()
            || !energy_residual_j.is_finite()
            || !volume_residual_m3.is_finite()
            || !entropy_change_j_k.is_finite()
            || stored_energy_change_j == 0.0
            || stored_energy_change_j.signum() != heat_j.signum()
            || energy_residual_j.abs() > tolerance_j
            || volume_residual_m3.abs() > tolerance_m3
            || !(self.mass_kg * solved.state().specific_internal_energy_j_kg()).is_finite()
            || solved.state().constitutive_law() != self.state.constitutive_law()
        {
            return Err(E::Invalid(
                "mixture parcel source energy or volume balance refused",
            ));
        }
        let report = MixtureEquilibriumHeatStep {
            supplied_heat_j: heat_j,
            stored_energy_change_j,
            energy_residual_j,
            volume_residual_m3,
            entropy_change_j_k,
            pressure_evaluations: solved.pressure_evaluations(),
            pt_evaluations: solved.pt_evaluations(),
            eos_evaluation_upper_bound: solved.eos_evaluation_upper_bound(),
        };
        poll(&mut cancelled)?;
        self.state = solved.into_state();
        Ok(report)
    }
}

fn admit_inventory(mass: f64, volume: f64) -> Result<f64, E> {
    let specific_volume = volume / mass;
    if !positive(mass)
        || !positive(volume)
        || !positive(specific_volume)
        || !positive(mass / volume)
    {
        return Err(E::Invalid(
            "positive representable mixture mass and volume required",
        ));
    }
    Ok(specific_volume)
}

fn source_specific_volume(state: &ResolvedThermochemicalPhaseState) -> f64 {
    state.molar_volume_m3_mol() / state.constitutive_law().molar_mass().value()
}
