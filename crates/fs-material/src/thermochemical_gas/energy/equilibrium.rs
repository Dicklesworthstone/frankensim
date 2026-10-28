//! Adiabatic constant-pressure equilibrium followed by an explicit chemistry freeze.
//!
//! The existing NASA-9 Gibbs-equilibrium owner supplies product composition and
//! temperature. Its complete species models must exactly match this material
//! bank, including source cards, molar masses and caloric corrections. Products
//! retain their original transport selections and pack identities.
//!
//! Formation enthalpy is already present in the source h(T); no separate heat
//! of reaction is added. The resulting parcel uses FROZEN composition for
//! storage, sound speed, heat capacity and downstream nozzle expansion. This
//! does not supply equilibrium nozzle chemistry, reaction rates, condensed
//! phases, ionization, a flame-stability verdict or a finite combustion chamber.

use super::parcel::ClosedGasParcel;
use super::{GasEnergyError, IsochoricGasParcel, updated_point};
use crate::thermochemical_gas::{ThermochemicalGasModel, coordinates, positive};
use fs_matdb::QueryPoint;
use fs_qty::{Amount, Energy, Temperature};
use fs_thermal::enthalpy::advection::species::close_generated_simplex;
use fs_thermochem::equilibrium::composition::{
    CompositionEquilibriumErrorV1, EquilibriumEnthalpySolveV1, IdealGasEquilibriumSolutionV1,
    IdealGasEquilibriumV1,
};
use fs_thermochem::{Composition, CompositionBasis, SpeciesId};

/// Explicit chamber chemistry limits and absolute final parcel-mass tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasEquilibriumChamberSolve {
    /// Source-temperature bracket, bounded work and total enthalpy tolerance [J].
    /// The same energy tolerance also bounds the final material parcel.
    pub equilibrium: EquilibriumEnthalpySolveV1,
    /// Absolute feed-to-parcel and equilibrium-to-parcel mass tolerance [kg].
    pub mass_tolerance_kg: f64,
}

/// Actual balances after composition canonicalization and material resolution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasEquilibriumChamberBalance {
    /// Sum of input amounts times their original source molar masses [kg].
    pub feed_mass_kg: f64,
    /// Final parcel mass minus feed mass [kg].
    pub mass_residual_kg: f64,
    /// Final parcel mass minus chemical-equilibrium product mass [kg].
    pub freeze_mass_residual_kg: f64,
    /// Actual final parcel mass times its resolved specific enthalpy [J].
    pub parcel_enthalpy_j: f64,
    /// Final parcel enthalpy minus the prescribed total feed enthalpy [J].
    pub enthalpy_residual_j: f64,
    /// Final parcel enthalpy minus the chemical-equilibrium enthalpy [J].
    pub freeze_enthalpy_residual_j: f64,
    /// Largest absolute mole-fraction change in exact-unit-sum roundoff closure.
    /// This is not a trace-species cutoff or a chemistry error estimate.
    pub maximum_mole_fraction_correction: f64,
}

/// Accepted equilibrium endpoint and its original-owner frozen gas parcel.
#[derive(Clone, Debug, PartialEq)]
pub struct GasEquilibriumChamber {
    equilibrium: IdealGasEquilibriumSolutionV1,
    parcel: ClosedGasParcel,
    balance: GasEquilibriumChamberBalance,
}
impl GasEquilibriumChamber {
    /// Reacting chamber composition, source receipts and nonlinear residuals.
    #[must_use]
    pub const fn equilibrium(&self) -> &IdealGasEquilibriumSolutionV1 {
        &self.equilibrium
    }
    /// Frozen products for existing storage, acoustic, thermal and nozzle APIs.
    #[must_use]
    pub const fn parcel(&self) -> &ClosedGasParcel {
        &self.parcel
    }
    /// Actual source-resolved mass, energy and canonicalization discrepancies.
    #[must_use]
    pub const fn balance(&self) -> GasEquilibriumChamberBalance {
        self.balance
    }
    /// Transfer ownership without re-resolving sources or adding reaction heat.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        ClosedGasParcel,
        IdealGasEquilibriumSolutionV1,
        GasEquilibriumChamberBalance,
    ) {
        (self.parcel, self.equilibrium, self.balance)
    }
}

/// Atomic chamber refusal; neither the material bank nor chemistry is mutated.
#[derive(Clone, Debug, PartialEq)]
pub enum GasEquilibriumChamberError {
    /// Invalid controls or unrepresentable derived physical quantities.
    Invalid(&'static str),
    /// Chemistry and transport bank do not retain the exact same species model.
    SourceMismatch(SpeciesId),
    /// A finite logarithmic abundance cannot be represented in the frozen model.
    /// Silently dropping this positive trace would change the admitted chemistry.
    UnrepresentableTrace(SpeciesId),
    /// The existing Gibbs-equilibrium owner refused the chamber state.
    Chemistry(CompositionEquilibriumErrorV1),
    /// The existing material/state owner refused the final frozen products.
    State(GasEnergyError),
    /// Actual final mass or energy did not meet the caller's absolute tolerance.
    Balance(GasEquilibriumChamberBalance),
    /// Cancellation before publishing any result.
    Cancelled,
}
impl core::fmt::Display for GasEquilibriumChamberError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermochemical equilibrium chamber: {self:?}")
    }
}
impl std::error::Error for GasEquilibriumChamberError {}
impl From<CompositionEquilibriumErrorV1> for GasEquilibriumChamberError {
    fn from(error: CompositionEquilibriumErrorV1) -> Self {
        match error {
            CompositionEquilibriumErrorV1::Cancelled => Self::Cancelled,
            other => Self::Chemistry(other),
        }
    }
}
impl From<GasEnergyError> for GasEquilibriumChamberError {
    fn from(error: GasEnergyError) -> Self {
        match error {
            GasEnergyError::Cancelled => Self::Cancelled,
            other => Self::State(other),
        }
    }
}
impl From<crate::MaterialError> for GasEquilibriumChamberError {
    fn from(error: crate::MaterialError) -> Self {
        Self::State(error.into())
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasEquilibriumChamberError> {
    if cancelled() {
        Err(GasEquilibriumChamberError::Cancelled)
    } else {
        Ok(())
    }
}

impl ThermochemicalGasModel {
    /// Adiabatic fixed-pressure equilibrium for feed species at one temperature.
    ///
    /// The context supplies the typed chamber pressure and any other material
    /// coordinates; its temperature is replaced by the solved chamber value.
    /// This model is a complete source bank: its catalog fractions are not the
    /// inlet composition. The feed declares actual species amounts in moles.
    /// Product volume follows from product mass and resolved density.
    ///
    /// For different inlet temperatures or prescribed stream enthalpy, use
    /// [Self::equilibrium_chamber_hp_enthalpy]. Formation energy is included
    /// once by the existing NASA-9 enthalpy owner.
    pub fn equilibrium_chamber_hp(
        &self,
        chemistry: &IdealGasEquilibriumV1,
        context: &QueryPoint,
        feed: &[(SpeciesId, Amount)],
        inlet_temperature: Temperature,
        controls: GasEquilibriumChamberSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEquilibriumChamber, GasEquilibriumChamberError> {
        self.admit_equilibrium_bank(chemistry, controls, &mut cancelled)?;
        let target = chemistry.feed_enthalpy(inlet_temperature, feed, &mut cancelled)?;
        self.equilibrium_chamber_hp_enthalpy(chemistry, context, feed, target, controls, cancelled)
    }

    /// Adiabatic fixed-pressure equilibrium for an explicit total enthalpy [J].
    ///
    /// Sum the independently sourced stream enthalpies in the SAME formation
    /// reference convention. Do not add reaction heat separately. The common
    /// material bank must contain every chemistry species, including potential
    /// products absent from the feed, with exact original NASA-9 model identity.
    ///
    /// Final products are frozen and resolved once through the original gas
    /// material owner. An underflowed positive trace, a missing transport source,
    /// an out-of-domain state or any failed balance refuses the whole result.
    pub fn equilibrium_chamber_hp_enthalpy(
        &self,
        chemistry: &IdealGasEquilibriumV1,
        context: &QueryPoint,
        feed: &[(SpeciesId, Amount)],
        target_enthalpy: Energy,
        controls: GasEquilibriumChamberSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEquilibriumChamber, GasEquilibriumChamberError> {
        self.admit_equilibrium_bank(chemistry, controls, &mut cancelled)?;
        let (_, pressure) = coordinates(context)?;
        let equilibrium = chemistry.equilibrate_hp_enthalpy(
            pressure,
            feed,
            target_enthalpy,
            controls.equilibrium,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;

        // Canonical source order makes this independent of caller feed ordering.
        // Chemistry has admitted each amount and rejected duplicate identities.
        let mut feed_mass = 0.0;
        let mut rounding = 0.0;
        for model in chemistry.models() {
            if let Some((_, amount)) = feed.iter().find(|(id, _)| id == model.species()) {
                let term = amount.value() * model.molar_mass().value();
                let adjusted = term - rounding;
                let next = feed_mass + adjusted;
                rounding = (next - feed_mass) - adjusted;
                feed_mass = next;
            }
        }
        if !positive(feed_mass) {
            return Err(GasEquilibriumChamberError::Invalid(
                "unrepresentable feed mass",
            ));
        }

        let mut fractions = Vec::with_capacity(equilibrium.species().len());
        for product in equilibrium.species() {
            if product.log_mole_fraction.is_some()
                && (!positive(product.mole_fraction) || !positive(product.amount.value()))
            {
                return Err(GasEquilibriumChamberError::UnrepresentableTrace(
                    product.species.clone(),
                ));
            }
            fractions.push(product.mole_fraction);
        }
        close_generated_simplex(&mut fractions).map_err(|error| {
            crate::thermochemical_gas::failure("equilibrium product simplex", error)
        })?;
        let correction = fractions
            .iter()
            .zip(equilibrium.species())
            .map(|(&fraction, product)| (fraction - product.mole_fraction).abs())
            .fold(0.0_f64, f64::max);
        let composition = Composition::new(CompositionBasis::MoleFraction, fractions)
            .map_err(|error| crate::thermochemical_gas::failure("equilibrium products", error))?;
        let model = self.with_composition(composition)?;
        let point = updated_point(
            context,
            equilibrium.temperature().value(),
            equilibrium.pressure().value(),
        )?;
        poll(&mut cancelled)?;
        let state = model.evaluate(&point)?;
        let volume = equilibrium.mass().value() / state.state().density;
        let mass = state.state().density * volume;
        if !positive(volume) || !positive(mass) {
            return Err(GasEquilibriumChamberError::Invalid(
                "unrepresentable chamber volume or mass",
            ));
        }
        let enthalpy = mass * state.specific_enthalpy_j_kg();
        let balance = GasEquilibriumChamberBalance {
            feed_mass_kg: feed_mass,
            mass_residual_kg: mass - feed_mass,
            freeze_mass_residual_kg: mass - equilibrium.mass().value(),
            parcel_enthalpy_j: enthalpy,
            enthalpy_residual_j: enthalpy - target_enthalpy.value(),
            freeze_enthalpy_residual_j: enthalpy - equilibrium.enthalpy().value(),
            maximum_mole_fraction_correction: correction,
        };
        if !enthalpy.is_finite()
            || balance.mass_residual_kg.abs() > controls.mass_tolerance_kg
            || balance.freeze_mass_residual_kg.abs() > controls.mass_tolerance_kg
            || !balance.enthalpy_residual_j.is_finite()
            || balance.enthalpy_residual_j.abs() > controls.equilibrium.enthalpy_tolerance_j
            || !balance.freeze_enthalpy_residual_j.is_finite()
            || balance.freeze_enthalpy_residual_j.abs() > controls.equilibrium.enthalpy_tolerance_j
        {
            return Err(GasEquilibriumChamberError::Balance(balance));
        }
        let parcel = ClosedGasParcel::from_isochoric(IsochoricGasParcel {
            model,
            point,
            state,
            mass_kg: mass,
            volume_m3: volume,
        });
        poll(&mut cancelled)?;
        Ok(GasEquilibriumChamber {
            equilibrium,
            parcel,
            balance,
        })
    }

    fn admit_equilibrium_bank(
        &self,
        chemistry: &IdealGasEquilibriumV1,
        controls: GasEquilibriumChamberSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), GasEquilibriumChamberError> {
        poll(cancelled)?;
        if !positive(controls.mass_tolerance_kg) {
            return Err(GasEquilibriumChamberError::Invalid(
                "positive finite chamber mass tolerance required",
            ));
        }
        let source = self.thermodynamics().components();
        let models = chemistry.models();
        if source.len() != models.len() {
            return Err(GasEquilibriumChamberError::Invalid(
                "material bank must cover the complete chemistry species set",
            ));
        }
        for (material, chemical) in source.iter().zip(models) {
            if material != chemical {
                return Err(GasEquilibriumChamberError::SourceMismatch(
                    chemical.species().clone(),
                ));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
