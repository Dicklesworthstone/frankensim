//! Closed mixture heating against a constant, quasistatic pressure load.
//!
//! The original phase-aware PH owner determines temperature, phase amounts and
//! volume. Heat is Delta H; the independently integrated Delta U and the actual
//! boundary work p Delta V must close the first law before publication. This
//! assumes mechanical equilibrium with a prescribed constant load. It does not
//! integrate piston inertia, heat-transfer rates or evaporation kinetics.

use super::{E, ThermochemicalPengRobinsonModel, poll, positive, source_specific_volume};
use crate::thermochemical_gas::nonideal::flash::equilibrium::{
    PhaseEnthalpyControls, ResolvedThermochemicalPhaseState,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::eos::peng_robinson::phase_screen::PengRobinsonPhaseScreenOptionsV1;

/// Accepted heat, expansion work and original source changes for a closed load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricMixtureHeatStep {
    /// Signed heat supplied to the mixture [J].
    pub supplied_heat_j: f64,
    /// Source-integrated enthalpy increment [J], with species datums cancelled.
    pub enthalpy_change_j: f64,
    /// Independently source-integrated internal-energy increment [J].
    pub stored_energy_change_j: f64,
    /// Work done BY the mixture against the constant pressure load [J].
    /// Positive expansion exports energy; compression imports it.
    pub boundary_work_j: f64,
    /// Delta U + p Delta V - Q on the actually published volume [J].
    pub energy_residual_j: f64,
    /// Source Delta H minus supplied heat [J].
    pub enthalpy_residual_j: f64,
    /// Published final volume minus initial volume [m3].
    pub volume_change_m3: f64,
    /// Source-integrated entropy change [J/K]. Heat alone does not specify the
    /// reservoir temperature, entropy transfer or process entropy production.
    pub entropy_change_j_k: f64,
    /// Actual original PT selections, including the reference and endpoints.
    pub pt_evaluations: usize,
    /// Original screening work plus split work ceilings, not observed EOS work.
    pub eos_evaluation_upper_bound: usize,
}

/// Conserved species inventory at one fixed prescribed external pressure.
///
/// Volume follows the original additive phase volumes. The original complete
/// model and phase evidence are retained; mass and composition cannot change.
/// Latent energy is part of the original thermodynamic state, never another
/// heat source. A clone preserves the full source and accepted state for retry.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricMixtureEquilibriumParcel {
    model: ThermochemicalPengRobinsonModel,
    state: ResolvedThermochemicalPhaseState,
    mass_kg: f64,
    pressure_pa: f64,
    volume_m3: f64,
}

impl IsobaricMixtureEquilibriumParcel {
    /// Resolve an initial temperature at the prescribed constant pressure.
    /// Initialization does not require subtraction of absolute energy datums.
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
            return Err(E::Invalid("positive finite isobaric mixture mass required"));
        }
        let state = model.evaluate_phase_equilibrium(
            temperature,
            pressure,
            phase_options,
            &mut cancelled,
        )?;
        let parcel = Self::checked(model, state, mass_kg)?;
        poll(&mut cancelled)?;
        Ok(parcel)
    }

    /// Resolve an initial source-reference enthalpy at the prescribed pressure.
    /// Use PT initialization plus heat increments when a large absolute datum
    /// cannot represent the requested caloric accuracy.
    pub fn from_enthalpy(
        model: ThermochemicalPengRobinsonModel,
        mass_kg: f64,
        pressure: Pressure,
        specific_enthalpy_j_kg: f64,
        controls: PhaseEnthalpyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, E> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !positive(mass_kg) {
            return Err(E::Invalid("positive finite isobaric mixture mass required"));
        }
        let solved = model.phase_equilibrium_pressure_enthalpy(
            pressure,
            specific_enthalpy_j_kg,
            controls,
            &mut cancelled,
        )?;
        let parcel = Self::checked(model, solved.into_state(), mass_kg)?;
        poll(&mut cancelled)?;
        Ok(parcel)
    }

    fn checked(
        model: ThermochemicalPengRobinsonModel,
        state: ResolvedThermochemicalPhaseState,
        mass_kg: f64,
    ) -> Result<Self, E> {
        let volume_m3 = admitted_volume(&state, mass_kg)?;
        Ok(Self {
            pressure_pa: state.pressure_pa(),
            model,
            state,
            mass_kg,
            volume_m3,
        })
    }

    /// Original immutable thermochemical and interaction sources.
    #[must_use]
    pub const fn model(&self) -> &ThermochemicalPengRobinsonModel {
        &self.model
    }
    /// Accepted phase state and original source/stability receipts.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalPhaseState {
        &self.state
    }
    /// Conserved total mass [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Prescribed constant mechanical pressure load [Pa].
    #[must_use]
    pub const fn pressure_pa(&self) -> f64 {
        self.pressure_pa
    }
    /// Actual volume from the accepted original EOS and phase amounts [m3].
    #[must_use]
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Accepted total internal energy on its original reference [J]. Use the
    /// increment report for small changes under large species energy datums.
    #[must_use]
    pub fn internal_energy_j(&self) -> f64 {
        self.mass_kg * self.state.specific_internal_energy_j_kg()
    }
    /// Accepted total enthalpy on its original reference [J].
    #[must_use]
    pub fn enthalpy_j(&self) -> f64 {
        self.mass_kg * self.state.specific_enthalpy_j_kg()
    }

    /// Add signed heat while the closed mixture expands against a constant load.
    ///
    /// The original PH solve can cross bubble/dew boundaries. Both source Delta H
    /// and Delta U cancel species reference datums before rounding. The reported
    /// work uses the volume that will actually be published, and the independent
    /// first-law check runs before either state or volume changes. Work/source
    /// failures, unrepresentable heat and final cancellation leave all fields
    /// unchanged. Zero heat preserves every accepted bit without a phase solve.
    pub fn add_heat(
        &mut self,
        heat_j: f64,
        controls: PhaseEnthalpyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricMixtureHeatStep, E> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !heat_j.is_finite() {
            return Err(E::Invalid("finite signed isobaric mixture heat required"));
        }
        if heat_j == 0.0 {
            let report = IsobaricMixtureHeatStep {
                supplied_heat_j: heat_j,
                enthalpy_change_j: 0.0,
                stored_energy_change_j: 0.0,
                boundary_work_j: 0.0,
                energy_residual_j: 0.0,
                enthalpy_residual_j: 0.0,
                volume_change_m3: 0.0,
                entropy_change_j_k: 0.0,
                pt_evaluations: 0,
                eos_evaluation_upper_bound: 0,
            };
            poll(&mut cancelled)?;
            return Ok(report);
        }
        let change_j_kg = heat_j / self.mass_kg;
        let tolerance_j = self.mass_kg * controls.enthalpy.enthalpy_tolerance_j_kg;
        if !change_j_kg.is_finite() || change_j_kg == 0.0 || !positive(tolerance_j) {
            return Err(E::Invalid(
                "isobaric mixture heat or tolerance is unrepresentable",
            ));
        }
        // A caller's coarse absolute tolerance must not swallow nonzero heat.
        let mut inner = controls;
        inner.enthalpy.enthalpy_tolerance_j_kg = inner
            .enthalpy
            .enthalpy_tolerance_j_kg
            .min(0.25 * change_j_kg.abs());
        inner.validate()?;
        let pressure = Pressure::new(self.pressure_pa);
        let solved = self.model.phase_equilibrium_pressure_enthalpy_increment(
            Temperature::new(self.state.temperature_k()),
            pressure,
            pressure,
            change_j_kg,
            inner,
            &mut cancelled,
        )?;
        let enthalpy_change_j = self.mass_kg
            * solved.accepted_increment_j_kg().ok_or(E::Invalid(
                "isobaric mixture source enthalpy increment is missing",
            ))?;
        let entropy_change_j_k = self.mass_kg
            * solved.entropy_change_j_kg_k().ok_or(E::Invalid(
                "isobaric mixture source entropy increment is missing",
            ))?;
        let volume_m3 = admitted_volume(solved.state(), self.mass_kg)?;
        let volume_change_m3 = volume_m3 - self.volume_m3;
        let boundary_work_j = self.pressure_pa * volume_change_m3;
        let stored_energy_change_j = self.mass_kg
            * self
                .model
                .phase_internal_energy_change(&self.state, solved.state())?;
        let enthalpy_residual_j = enthalpy_change_j - heat_j;
        let energy_residual_j = (stored_energy_change_j - heat_j) + boundary_work_j;
        if [
            enthalpy_change_j,
            entropy_change_j_k,
            volume_change_m3,
            boundary_work_j,
            stored_energy_change_j,
            enthalpy_residual_j,
            energy_residual_j,
        ]
        .into_iter()
        .any(|value| !value.is_finite())
            || enthalpy_change_j == 0.0
            || enthalpy_change_j.signum() != heat_j.signum()
            || enthalpy_residual_j.abs() > tolerance_j
            || energy_residual_j.abs() > tolerance_j
            || solved.state().pressure_pa().to_bits() != self.pressure_pa.to_bits()
        {
            return Err(E::Invalid(
                "isobaric mixture original energy/work balance failed",
            ));
        }
        let report = IsobaricMixtureHeatStep {
            supplied_heat_j: heat_j,
            enthalpy_change_j,
            stored_energy_change_j,
            boundary_work_j,
            energy_residual_j,
            enthalpy_residual_j,
            volume_change_m3,
            entropy_change_j_k,
            pt_evaluations: solved.pt_evaluations(),
            eos_evaluation_upper_bound: solved.eos_evaluation_upper_bound(),
        };
        poll(&mut cancelled)?;
        self.state = solved.into_state();
        self.volume_m3 = volume_m3;
        Ok(report)
    }
}

fn admitted_volume(state: &ResolvedThermochemicalPhaseState, mass: f64) -> Result<f64, E> {
    let volume = mass * source_specific_volume(state);
    if !positive(mass)
        || !positive(volume)
        || !positive(mass / volume)
        || !(mass * state.specific_internal_energy_j_kg()).is_finite()
        || !(mass * state.specific_enthalpy_j_kg()).is_finite()
    {
        Err(E::Invalid("isobaric mixture inventory is unrepresentable"))
    } else {
        Ok(volume)
    }
}
