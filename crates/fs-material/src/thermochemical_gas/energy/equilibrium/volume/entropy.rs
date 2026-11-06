//! Entropy-constrained closed-gas volume changes with source-preserving products.
//!
//! The original reacting SV owner supplies temperature, pressure and composition.
//! The common material boundary then retains exact feed mass and target volume,
//! original calorics, transport cards and source receipts. Entropy is checked
//! again after composition rounding and the mass-preserving pressure correction.
//!
//! This is a reversible equilibrium endpoint when starting from an equilibrated
//! state. It is not a finite-rate piston, ignition or shock model. Returned
//! material products are explicitly frozen for subsequent transport and storage;
//! reacting response remains available from the accompanying equilibrium state.

use super::{
    FreezeVolumeControls, GasEquilibriumChamberError, GasIsochoricEquilibriumBalance,
    GasIsochoricEquilibriumChamber, GasIsochoricEquilibriumError,
    IdealGasIsochoricEquilibriumV1, ThermochemicalGasModel, poll, positive,
};
use crate::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_matdb::QueryPoint;
use fs_qty::{Amount, Volume};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::composition::{
    IdealGasEquilibriumV1,
    volume::{EquilibriumVolumeErrorV1, entropy::EquilibriumVolumeEntropySolveV1},
};

/// Independent chemical, material-mass and composition-freezing tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasIsentropicVolumeSolve {
    /// Original chemical-equilibrium, volume and total-entropy solve controls.
    pub equilibrium: EquilibriumVolumeEntropySolveV1,
    /// Maximum absolute feed/product/density mass discrepancy [kg].
    pub mass_tolerance_kg: f64,
    /// Maximum absolute internal-energy change from freezing the products [J].
    /// Energy is an output of SV, not a separately prescribed conserved target.
    pub freeze_energy_tolerance_j: f64,
}

/// Actual material-state discrepancies at the accepted SV endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasIsentropicVolumeBalance {
    /// Common source-preserving mass, pressure and energy-freezing balance.
    /// Its energy residual is relative to the native SV endpoint's output U.
    pub freeze: GasIsochoricEquilibriumBalance,
    /// Entropy from the final parcel's original species-plus-mixing law [J/K].
    pub parcel_entropy_j_k: f64,
    /// Final material entropy minus the requested total entropy [J/K].
    pub entropy_residual_j_k: f64,
    /// Final material entropy minus the native reacting endpoint entropy [J/K].
    pub freeze_entropy_residual_j_k: f64,
}

/// Reacting SV endpoint and its source-preserving frozen material products.
#[derive(Clone, Debug, PartialEq)]
pub struct GasIsentropicVolumeChamber {
    chamber: GasIsochoricEquilibriumChamber,
    balance: GasIsentropicVolumeBalance,
}
impl GasIsentropicVolumeChamber {
    /// Native composition, entropy, reacting response and convergence residuals.
    #[must_use]
    pub const fn equilibrium(&self) -> &IdealGasIsochoricEquilibriumV1 {
        self.chamber.equilibrium()
    }
    /// Frozen products for existing thermal, acoustic and mechanical consumers.
    #[must_use]
    pub const fn parcel(&self) -> &ClosedGasParcel {
        self.chamber.parcel()
    }
    /// Actual balances after source resolution and exact-mass EOS closure.
    #[must_use]
    pub const fn balance(&self) -> GasIsentropicVolumeBalance {
        self.balance
    }
    /// Transfer accepted products and evidence without another equilibrium solve.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        ClosedGasParcel,
        IdealGasIsochoricEquilibriumV1,
        GasIsentropicVolumeBalance,
    ) {
        let (parcel, equilibrium, _) = self.chamber.into_parts();
        (parcel, equilibrium, self.balance)
    }
}

/// Atomic refusal; no input material bank or existing parcel is modified.
#[derive(Clone, Debug, PartialEq)]
pub enum GasIsentropicVolumeError {
    /// Original chemistry, source, trace or mass/energy-freezing refusal.
    Volume(GasIsochoricEquilibriumError),
    /// Actual final entropy or entropy-freezing discrepancy exceeds its limit.
    EntropyBalance(GasIsentropicVolumeBalance),
    /// Invalid extra material-freezing controls.
    Invalid(&'static str),
    /// Cancellation before publishing any result.
    Cancelled,
}
impl core::fmt::Display for GasIsentropicVolumeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "isentropic volume equilibrium chamber: {self:?}")
    }
}
impl std::error::Error for GasIsentropicVolumeError {}
impl From<GasIsochoricEquilibriumError> for GasIsentropicVolumeError {
    fn from(error: GasIsochoricEquilibriumError) -> Self {
        match error {
            GasIsochoricEquilibriumError::Cancelled => Self::Cancelled,
            other => Self::Volume(other),
        }
    }
}
impl From<GasEquilibriumChamberError> for GasIsentropicVolumeError {
    fn from(error: GasEquilibriumChamberError) -> Self {
        GasIsochoricEquilibriumError::from(error).into()
    }
}
impl From<EquilibriumVolumeErrorV1> for GasIsentropicVolumeError {
    fn from(error: EquilibriumVolumeErrorV1) -> Self {
        GasIsochoricEquilibriumError::from(error).into()
    }
}

impl ThermochemicalGasModel {
    /// Resolve a closed gas at a prescribed total entropy, volume and inventory.
    ///
    /// Every chemical temperature trial uses reacting composition and original
    /// variable heat capacity. The material bank must cover the complete neutral
    /// species set with exact original caloric identities; active charged species
    /// still require the separate ionized transport owner. The context supplies
    /// other material axes; its temperature and pressure are replaced.
    ///
    /// Total internal energy changes during reversible compression/expansion.
    /// For an adiabatic path between equilibrated endpoints, work on the gas is
    /// the difference of the returned total internal energies. No extra heat of
    /// reaction is added. Source escape, trace underflow, failed balances or
    /// cancellation returns no partially resolved product parcel.
    #[allow(clippy::too_many_arguments)]
    pub fn equilibrium_chamber_sv(
        &self,
        chemistry: &IdealGasEquilibriumV1,
        context: &QueryPoint,
        feed: &[(SpeciesId, Amount)],
        target_entropy_j_k: f64,
        volume: Volume,
        controls: GasIsentropicVolumeSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasIsentropicVolumeChamber, GasIsentropicVolumeError> {
        poll(&mut cancelled)?;
        if !positive(controls.freeze_energy_tolerance_j) {
            return Err(GasIsentropicVolumeError::Invalid(
                "positive finite SV product energy-freezing tolerance required",
            ));
        }
        self.admit_equilibrium_bank(chemistry, controls.mass_tolerance_kg, &mut cancelled)?;
        let equilibrium = chemistry.equilibrate_sv(
            volume,
            feed,
            target_entropy_j_k,
            controls.equilibrium,
            &mut cancelled,
        )?;
        let chamber = self.freeze_volume_equilibrium(
            chemistry,
            context,
            feed,
            equilibrium,
            None,
            FreezeVolumeControls {
                mass_tolerance_kg: controls.mass_tolerance_kg,
                energy_tolerance_j: controls.freeze_energy_tolerance_j,
                relative_volume_tolerance: controls.equilibrium.volume.relative_volume_tolerance,
            },
            &mut cancelled,
        )?;
        let properties = chamber.parcel().state().thermodynamics().properties();
        let entropy = (chamber.parcel().mass_kg() / properties.molar_mass().value())
            * properties.s().value();
        let balance = GasIsentropicVolumeBalance {
            freeze: chamber.balance(),
            parcel_entropy_j_k: entropy,
            entropy_residual_j_k: entropy - target_entropy_j_k,
            freeze_entropy_residual_j_k: entropy
                - chamber.equilibrium().thermodynamic().entropy_j_k(),
        };
        let tolerance = controls.equilibrium.entropy_tolerance_j_k;
        if !entropy.is_finite()
            || !balance.entropy_residual_j_k.is_finite()
            || !balance.freeze_entropy_residual_j_k.is_finite()
            || balance.entropy_residual_j_k.abs() > tolerance
            || balance.freeze_entropy_residual_j_k.abs() > tolerance
        {
            return Err(GasIsentropicVolumeError::EntropyBalance(balance));
        }
        poll(&mut cancelled)?;
        Ok(GasIsentropicVolumeChamber { chamber, balance })
    }
}

#[cfg(test)]
mod tests;
