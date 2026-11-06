//! A moving closed parcel on one prepared reacting equilibrium isentrope.
//!
//! The original species/element inventory and reference entropy never change.
//! Density queries reuse the native Gibbs workspace; only accepted endpoints
//! resolve source-preserving material products. Relative source calorics supply
//! work without subtracting large absolute formation energies.
//!
//! This is the reversible, infinitely fast neutral ideal-gas chemistry limit.
//! Frozen products remain available for transport, while pressure stiffness is
//! explicitly the zero-frequency reacting response. There is no finite-rate
//! piston, irreversible-work, shock or frequency-dependent relaxation model.
use super::{
    FreezeVolumeControls, GasIsochoricEquilibriumChamber, GasIsochoricEquilibriumError,
    ThermochemicalGasModel,
};
use crate::thermochemical_gas::energy::{
    entropy::GasEntropyError,
    parcel::mechanics::{GasMechanicalLoad, GasPressureResponse, MAX_VOLUME_COORDINATES},
};
use fs_matdb::QueryPoint;
use fs_qty::{Amount, Pressure, Temperature, Volume};
use fs_thermochem::equilibrium::composition::{
    CompositionEquilibriumErrorV1, EquilibriumEntropySolveV1, EquilibriumExpansionErrorV1,
    IdealGasEquilibriumV1,
    isentrope::{
        EquilibriumIsentropeSampleV1, EquilibriumIsentropeV1,
        density::{EquilibriumDensityErrorV1, EquilibriumDensitySolveV1},
    },
    volume::IdealGasIsochoricEquilibriumV1,
};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

/// Fixed source domains, bounded query work and whole-path conservation limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumGasParcelSolve {
    /// Original composition/temperature limits and extensive entropy tolerance.
    pub thermochemistry: EquilibriumEntropySolveV1,
    /// Admitted pressures and density/volume accuracy at every endpoint.
    pub density: EquilibriumDensitySolveV1,
    /// Feed/product/material mass discrepancy limit [kg].
    pub mass_tolerance_kg: f64,
    /// Absolute energy-freezing, step balance and reference balance limit [J].
    /// This is one allowance at any endpoint, never multiplied by step count.
    pub energy_tolerance_j: f64,
}

/// Actual change between two accepted reversible volume endpoints.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumGasVolumeChange {
    /// Previous exact parcel volume [m3].
    pub previous_volume_m3: f64,
    /// New exact parcel volume [m3].
    pub volume_m3: f64,
    /// Actual frozen-material internal-energy difference [J].
    pub stored_energy_change_j: f64,
    /// Reversible work OUT of the gas, from source-relative calorics [J].
    pub work_by_gas_j: f64,
    /// Actual stored energy change plus work by gas [J].
    pub energy_residual_j: f64,
    /// Work OUT since the immutable birth state [J], without a running sum.
    pub work_by_gas_since_reference_j: f64,
    /// Actual storage change since birth plus work since birth [J].
    pub reference_energy_residual_j: f64,
    /// Material entropy change over this update [J/K].
    pub entropy_change_j_k: f64,
    /// Material entropy change since birth [J/K], on the original reference.
    pub entropy_change_since_reference_j_k: f64,
    /// Pressure queries in this update; zero for an exact volume no-op.
    pub pressure_queries: usize,
    /// Accumulated composition iterations, including discarded pressure trials.
    pub composition_iterations: usize,
    /// Accumulated temperature iterations, including discarded pressure trials.
    pub temperature_iterations: usize,
}

/// No failed query, source resolution or balance check publishes a new parcel.
#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumGasParcelError {
    /// Original chemistry, temperature-domain or entropy refusal.
    Thermochemistry(EquilibriumExpansionErrorV1),
    /// Original density/pressure/volume query refusal.
    Density(EquilibriumDensityErrorV1),
    /// Original material source, mass or caloric-freezing refusal.
    Material(GasIsochoricEquilibriumError),
    /// Invalid or unrepresentable state/control.
    Invalid(&'static str),
    /// Actual step or whole-path energy/entropy failed the declared allowance.
    Balance(EquilibriumGasVolumeChange),
    /// Cancellation before publishing accepted physical state.
    Cancelled,
}
impl core::fmt::Display for EquilibriumGasParcelError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reversible reacting parcel: {self:?}")
    }
}
impl std::error::Error for EquilibriumGasParcelError {}
impl From<EquilibriumExpansionErrorV1> for EquilibriumGasParcelError {
    fn from(error: EquilibriumExpansionErrorV1) -> Self {
        match error {
            EquilibriumExpansionErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled) => {
                Self::Cancelled
            }
            other => Self::Thermochemistry(other),
        }
    }
}
impl From<EquilibriumDensityErrorV1> for EquilibriumGasParcelError {
    fn from(error: EquilibriumDensityErrorV1) -> Self {
        match error {
            EquilibriumDensityErrorV1::Thermochemistry(error) => error.into(),
            other => Self::Density(other),
        }
    }
}
impl From<GasIsochoricEquilibriumError> for EquilibriumGasParcelError {
    fn from(error: GasIsochoricEquilibriumError) -> Self {
        match error {
            GasIsochoricEquilibriumError::Cancelled => Self::Cancelled,
            other => Self::Material(other),
        }
    }
}
type Error = EquilibriumGasParcelError;
fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

/// One accepted physical chamber and one reusable numerical workspace.
///
/// The original bank and chemistry are borrowed immutably. Solver warm starts
/// may change on a refused trial, but the physical state, work and entropy
/// reference remain unchanged. A mechanical view borrows the accepted chamber,
/// so it cannot survive a mutable volume update.
pub struct EquilibriumGasParcel<'a> {
    bank: &'a ThermochemicalGasModel,
    chemistry: &'a IdealGasEquilibriumV1,
    feed: Vec<(SpeciesId, Amount)>,
    context: QueryPoint,
    path: EquilibriumIsentropeV1<'a>,
    controls: EquilibriumGasParcelSolve,
    chamber: GasIsochoricEquilibriumChamber,
    sample: EquilibriumIsentropeSampleV1,
    reference_temperature_k: f64,
    reference_amount_mol: f64,
    reference_energy_j: f64,
    reference_freeze_entropy_j_k: f64,
    energy_from_reference_j: f64,
    entropy_from_reference_j_k: f64,
}
impl core::fmt::Debug for EquilibriumGasParcel<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("EquilibriumGasParcel")
            .field("chamber", &self.chamber)
            .field("controls", &self.controls)
            .field("energy_from_reference_j", &self.energy_from_reference_j)
            .field(
                "entropy_from_reference_j_k",
                &self.entropy_from_reference_j_k,
            )
            .finish_non_exhaustive()
    }
}

impl<'a> EquilibriumGasParcel<'a> {
    /// Establish a reacting equilibrium at explicit T,p and feed inventory.
    ///
    /// Initial volume follows conserved feed mass and the original equilibrium
    /// density. The entire pressure/temperature domain is caller supplied.
    /// The feed specifies elemental inventory; the bank's catalog fractions
    /// do not freeze the initial or subsequent chemistry.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bank: &'a ThermochemicalGasModel,
        chemistry: &'a IdealGasEquilibriumV1,
        context: &QueryPoint,
        feed: &[(SpeciesId, Amount)],
        temperature: Temperature,
        pressure: Pressure,
        controls: EquilibriumGasParcelSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        checkpoint(&mut cancelled)?;
        if !controls.energy_tolerance_j.is_finite() || controls.energy_tolerance_j <= 0.0 {
            return Err(Error::Invalid(
                "positive finite reacting-parcel energy allowance required",
            ));
        }
        bank.admit_equilibrium_bank(chemistry, controls.mass_tolerance_kg, &mut cancelled)
            .map_err(GasIsochoricEquilibriumError::from)?;
        let mut path = EquilibriumIsentropeV1::new(
            chemistry,
            temperature,
            pressure,
            feed,
            controls.thermochemistry,
            &mut cancelled,
        )?;
        let reference_temperature_k = path.reference_state().equilibrium().temperature().value();
        let reference_amount_mol = path.reference_state().equilibrium().total_amount().value();
        let mass = path.feed_mass_kg();
        let volume = Volume::new(mass / path.reference_state().response().density_kg_m3);
        let (state, query) = path
            .sample_volume(volume, controls.density, &mut cancelled)?
            .into_parts();
        let (chamber, freeze_entropy) = freeze(
            bank,
            chemistry,
            context,
            feed,
            state,
            mass,
            controls,
            &mut cancelled,
        )?;
        let reference_energy_j = chamber.balance().parcel_internal_energy_j;
        checkpoint(&mut cancelled)?;
        Ok(Self {
            bank,
            chemistry,
            feed: feed.to_vec(),
            context: context.clone(),
            path,
            controls,
            chamber,
            sample: query.sample,
            reference_temperature_k,
            reference_amount_mol,
            reference_energy_j,
            reference_freeze_entropy_j_k: freeze_entropy,
            energy_from_reference_j: 0.0,
            entropy_from_reference_j_k: 0.0,
        })
    }

    /// Accepted source-bearing reacting endpoint and its frozen transport view.
    pub const fn chamber(&self) -> &GasIsochoricEquilibriumChamber {
        &self.chamber
    }
    /// Accepted reacting derivatives and direct source-relative increments.
    pub const fn sample(&self) -> EquilibriumIsentropeSampleV1 {
        self.sample
    }
    /// Unchanging original elemental/species feed inventory.
    pub fn feed(&self) -> &[(SpeciesId, Amount)] {
        &self.feed
    }
    /// Exact conserved feed mass [kg].
    pub fn mass_kg(&self) -> f64 {
        self.chamber.parcel().mass_kg()
    }
    /// Exact accepted chamber volume [m3].
    pub fn volume_m3(&self) -> f64 {
        self.chamber.parcel().volume_m3()
    }
    /// Fixed solver/source/work controls.
    pub const fn controls(&self) -> EquilibriumGasParcelSolve {
        self.controls
    }
    /// Reversible work OUT since birth, computed from one immutable reference.
    pub fn work_by_gas_since_reference_j(&self) -> f64 {
        -self.energy_from_reference_j
    }

    /// Re-equilibrate at a new volume and publish only after all physical checks.
    ///
    /// Temperature and composition both evolve. Work is delta H - delta(NRT),
    /// using the original path's direct caloric differences. It does not
    /// subtract large absolute formation energies or integrate a guessed gamma.
    /// Energy/entropy checks cover both this step AND the complete path.
    pub fn set_volume_isentropic(
        &mut self,
        volume: Volume,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumGasVolumeChange, Error> {
        checkpoint(&mut cancelled)?;
        let old_volume = self.volume_m3();
        let old_energy = self.chamber.balance().parcel_internal_energy_j;
        if volume.value() == old_volume {
            return Ok(EquilibriumGasVolumeChange {
                previous_volume_m3: old_volume,
                volume_m3: old_volume,
                stored_energy_change_j: 0.0,
                work_by_gas_j: 0.0,
                energy_residual_j: 0.0,
                work_by_gas_since_reference_j: -self.energy_from_reference_j,
                reference_energy_residual_j: (old_energy - self.reference_energy_j)
                    - self.energy_from_reference_j,
                entropy_change_j_k: 0.0,
                entropy_change_since_reference_j_k: self.entropy_from_reference_j_k,
                pressure_queries: 0,
                composition_iterations: 0,
                temperature_iterations: 0,
            });
        }
        let (state, query) = self
            .path
            .sample_volume(volume, self.controls.density, &mut cancelled)?
            .into_parts();
        let chemical = state.equilibrium();
        // Difference of N*T formed before multiplying R; reference formation
        // offsets never enter this pressure-volume increment.
        let n = chemical.total_amount().value();
        let delta_nrt = R
            * (n * (chemical.temperature().value() - self.reference_temperature_k)
                + self.reference_temperature_k * (n - self.reference_amount_mol));
        let energy_from_reference =
            self.mass_kg() * query.sample.specific_enthalpy_change_j_kg - delta_nrt;
        let (chamber, freeze_entropy) = freeze(
            self.bank,
            self.chemistry,
            &self.context,
            &self.feed,
            state,
            self.mass_kg(),
            self.controls,
            &mut cancelled,
        )?;
        let entropy_from_reference =
            query.sample.entropy_change_j_k + (freeze_entropy - self.reference_freeze_entropy_j_k);
        let energy = chamber.balance().parcel_internal_energy_j;
        let stored = energy - old_energy;
        let work = self.energy_from_reference_j - energy_from_reference;
        let report = EquilibriumGasVolumeChange {
            previous_volume_m3: old_volume,
            volume_m3: volume.value(),
            stored_energy_change_j: stored,
            work_by_gas_j: work,
            energy_residual_j: stored + work,
            work_by_gas_since_reference_j: -energy_from_reference,
            reference_energy_residual_j: (energy - self.reference_energy_j) - energy_from_reference,
            entropy_change_j_k: entropy_from_reference - self.entropy_from_reference_j_k,
            entropy_change_since_reference_j_k: entropy_from_reference,
            pressure_queries: query.pressure_queries,
            composition_iterations: query.composition_iterations,
            temperature_iterations: query.temperature_iterations,
        };
        if [
            stored,
            work,
            report.energy_residual_j,
            energy_from_reference,
            report.reference_energy_residual_j,
            report.entropy_change_j_k,
            entropy_from_reference,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || report.energy_residual_j.abs() > self.controls.energy_tolerance_j
            || report.reference_energy_residual_j.abs() > self.controls.energy_tolerance_j
            || report.entropy_change_j_k.abs() > self.controls.thermochemistry.entropy_tolerance_j_k
            || entropy_from_reference.abs() > self.controls.thermochemistry.entropy_tolerance_j_k
            || chamber.parcel().mass_kg() != self.mass_kg()
        {
            return Err(Error::Balance(report));
        }
        checkpoint(&mut cancelled)?;
        self.chamber = chamber;
        self.sample = query.sample;
        self.energy_from_reference_j = energy_from_reference;
        self.entropy_from_reference_j_k = entropy_from_reference;
        Ok(report)
    }

    /// Zero-frequency reacting pressure stiffness on the accepted material state.
    ///
    /// Gamma_s comes from the original differentiated chemical KKT equations.
    /// Multiplying it by the exact-mass material pressure retains the admitted
    /// pressure-freezing correction; that discrepancy remains in chamber.balance().
    pub fn pressure_response(&self) -> Result<GasPressureResponse, GasEntropyError> {
        let pressure = self.chamber.parcel().state().state().pressure;
        let bulk = pressure * self.sample.response.isentropic_pressure_exponent;
        let derivative = -bulk / self.volume_m3();
        if !pressure.is_finite()
            || pressure <= 0.0
            || !bulk.is_finite()
            || bulk <= 0.0
            || !derivative.is_finite()
            || derivative >= 0.0
        {
            return Err(GasEntropyError::Invalid(
                "unrepresentable reacting mechanical response",
            ));
        }
        Ok(GasPressureResponse {
            pressure_pa: pressure,
            bulk_modulus_pa: bulk,
            pressure_volume_derivative_pa_m3: derivative,
        })
    }

    /// Borrow work-conjugate pressure forces and matrix-free geometric tangents.
    pub fn mechanical_load<'b>(
        &'b self,
        volume_m3: f64,
        gradient: &'b [f64],
        exterior_pressure_pa: f64,
    ) -> Result<GasMechanicalLoad<'b, GasIsochoricEquilibriumChamber>, GasEntropyError> {
        self.mechanical_load_with_limit(
            volume_m3,
            gradient,
            exterior_pressure_pa,
            MAX_VOLUME_COORDINATES,
            || false,
        )
    }

    /// The same borrowed reacting response with caller-owned coordinate limits.
    pub fn mechanical_load_with_limit<'b>(
        &'b self,
        volume_m3: f64,
        gradient: &'b [f64],
        exterior_pressure_pa: f64,
        maximum_coordinates: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasMechanicalLoad<'b, GasIsochoricEquilibriumChamber>, GasEntropyError> {
        GasMechanicalLoad::new_with_limit(
            &self.chamber,
            self.volume_m3(),
            volume_m3,
            gradient,
            exterior_pressure_pa,
            self.pressure_response()?,
            maximum_coordinates,
            cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn freeze(
    bank: &ThermochemicalGasModel,
    chemistry: &IdealGasEquilibriumV1,
    context: &QueryPoint,
    feed: &[(SpeciesId, Amount)],
    state: IdealGasIsochoricEquilibriumV1,
    mass: f64,
    controls: EquilibriumGasParcelSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(GasIsochoricEquilibriumChamber, f64), Error> {
    let chamber = bank.freeze_volume_equilibrium(
        chemistry,
        context,
        feed,
        state,
        Some(mass),
        FreezeVolumeControls {
            mass_tolerance_kg: controls.mass_tolerance_kg,
            energy_tolerance_j: controls.energy_tolerance_j,
            relative_volume_tolerance: controls.density.relative_density_tolerance,
        },
        &mut *cancelled,
    )?;
    let properties = chamber.parcel().state().thermodynamics().properties();
    let entropy = (mass / properties.molar_mass().value()) * properties.s().value();
    let residual = entropy - chamber.equilibrium().thermodynamic().entropy_j_k();
    if !entropy.is_finite()
        || !residual.is_finite()
        || residual.abs() > controls.thermochemistry.entropy_tolerance_j_k
    {
        return Err(Error::Invalid(
            "reacting material entropy freezing exceeded its allowance",
        ));
    }
    checkpoint(cancelled)?;
    Ok((chamber, residual))
}

#[cfg(test)]
mod tests;
