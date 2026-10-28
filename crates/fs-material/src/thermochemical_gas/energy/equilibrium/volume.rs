//! Sealed, rigid ideal-gas equilibrium with one conserved internal energy.
//!
//! Chemistry changes molecular count, pressure and temperature while retaining
//! elemental inventory, mass and volume. The original NASA-9 formation energies
//! are included once. Products then use the existing frozen material owner for
//! transport, storage and acoustics; the reacting equilibrium response remains
//! separately available. No reaction rate, ignition or flame model is implied.

use super::{GasEquilibriumChamberError, ThermochemicalGasModel, equilibrium_feed_mass, poll};
use crate::thermochemical_gas::energy::{
    IsochoricGasParcel, parcel::ClosedGasParcel, updated_point,
};
use crate::thermochemical_gas::{failure, positive};
use fs_matdb::QueryPoint;
use fs_qty::{Amount, Energy, Temperature, Volume};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::composition::{
    CompositionEquilibriumErrorV1, EquilibriumExpansionErrorV1, IdealGasEquilibriumV1,
    volume::{
        EquilibriumInternalEnergySolveV1, EquilibriumVolumeErrorV1, IdealGasIsochoricEquilibriumV1,
    },
};

/// Absolute conservation limits for a source-resolved sealed chamber.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasIsochoricEquilibriumSolve {
    /// Chemical/volume closure and total internal-energy tolerance [J].
    pub equilibrium: EquilibriumInternalEnergySolveV1,
    /// Absolute feed, chemistry and final density/mass discrepancy limit [kg].
    pub mass_tolerance_kg: f64,
}

/// Actual balances after freezing composition and resolving all material data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasIsochoricEquilibriumBalance {
    /// Original feed mass [kg], reconstructed from the original species data.
    pub feed_mass_kg: f64,
    /// Preserved parcel mass minus feed mass [kg].
    pub mass_residual_kg: f64,
    /// Chemical-equilibrium product mass minus preserved parcel mass [kg].
    pub chemical_mass_residual_kg: f64,
    /// Resolved density times fixed volume minus preserved mass [kg].
    pub density_mass_residual_kg: f64,
    /// Actual stored internal energy of the resulting frozen parcel [J].
    pub parcel_internal_energy_j: f64,
    /// Actual stored internal energy minus the prescribed target [J].
    pub energy_residual_j: f64,
    /// Frozen-parcel internal energy minus the chemical-equilibrium value [J].
    pub freeze_energy_residual_j: f64,
    /// Relative pressure adjustment needed to retain the exact mass and volume.
    /// Bounded by the admitted volume/mass errors; not an extra flow process.
    pub relative_pressure_correction: f64,
    /// Largest mole-fraction change in the existing exact-unit-sum closure.
    pub maximum_mole_fraction_correction: f64,
}

/// Original reacting endpoint and a common source-preserving frozen gas parcel.
#[derive(Clone, Debug, PartialEq)]
pub struct GasIsochoricEquilibriumChamber {
    equilibrium: IdealGasIsochoricEquilibriumV1,
    parcel: ClosedGasParcel,
    balance: GasIsochoricEquilibriumBalance,
}

impl GasIsochoricEquilibriumChamber {
    /// Original equilibrium species, sources, residuals and reacting response.
    #[must_use]
    pub const fn equilibrium(&self) -> &IdealGasIsochoricEquilibriumV1 {
        &self.equilibrium
    }
    /// Frozen gas for existing thermal, mechanical, transport and acoustic APIs.
    #[must_use]
    pub const fn parcel(&self) -> &ClosedGasParcel {
        &self.parcel
    }
    /// Actual conservation and freezing discrepancies.
    #[must_use]
    pub const fn balance(&self) -> GasIsochoricEquilibriumBalance {
        self.balance
    }
    /// Transfer the accepted state without resampling material sources.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        ClosedGasParcel,
        IdealGasIsochoricEquilibriumV1,
        GasIsochoricEquilibriumBalance,
    ) {
        (self.parcel, self.equilibrium, self.balance)
    }
}

/// Atomic refusal: no material bank or existing parcel is modified.
#[derive(Clone, Debug, PartialEq)]
pub enum GasIsochoricEquilibriumError {
    /// Original chemistry, source-path or bounded solver refusal.
    Chemistry(EquilibriumVolumeErrorV1),
    /// Shared exact-source, trace or material-state admission refused.
    Material(GasEquilibriumChamberError),
    /// Actual final conservation or freezing errors exceeded the caller's limits.
    Balance(GasIsochoricEquilibriumBalance),
    /// Invalid or unrepresentable extensive state.
    Invalid(&'static str),
    /// Cancellation before publishing any result or parcel update.
    Cancelled,
}
impl core::fmt::Display for GasIsochoricEquilibriumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "isochoric equilibrium chamber: {self:?}")
    }
}
impl std::error::Error for GasIsochoricEquilibriumError {}
impl From<GasEquilibriumChamberError> for GasIsochoricEquilibriumError {
    fn from(error: GasEquilibriumChamberError) -> Self {
        match error {
            GasEquilibriumChamberError::Cancelled => Self::Cancelled,
            other => Self::Material(other),
        }
    }
}
impl From<EquilibriumVolumeErrorV1> for GasIsochoricEquilibriumError {
    fn from(error: EquilibriumVolumeErrorV1) -> Self {
        match error {
            EquilibriumVolumeErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled)
            | EquilibriumVolumeErrorV1::Response(EquilibriumExpansionErrorV1::Equilibrium(
                CompositionEquilibriumErrorV1::Cancelled,
            )) => Self::Cancelled,
            other => Self::Chemistry(other),
        }
    }
}

impl ThermochemicalGasModel {
    /// Adiabatic fixed-volume equilibrium for a feed at one initial temperature.
    ///
    /// The typed context supplies material axes; its temperature and pressure
    /// are replaced by the solved state. The complete bank must include possible
    /// products even when absent from the feed. Pressure follows molecular count
    /// and temperature; constant-pressure enthalpy closure is not used here.
    #[allow(clippy::too_many_arguments)]
    pub fn equilibrium_chamber_uv(
        &self,
        chemistry: &IdealGasEquilibriumV1,
        context: &QueryPoint,
        feed: &[(SpeciesId, Amount)],
        inlet_temperature: Temperature,
        volume: Volume,
        controls: GasIsochoricEquilibriumSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasIsochoricEquilibriumChamber, GasIsochoricEquilibriumError> {
        self.admit_equilibrium_bank(chemistry, controls.mass_tolerance_kg, &mut cancelled)?;
        let target = chemistry
            .feed_internal_energy(inlet_temperature, feed, &mut cancelled)
            .map_err(EquilibriumVolumeErrorV1::from)?;
        self.equilibrium_chamber_uv_energy(
            chemistry, context, feed, target, volume, controls, cancelled,
        )
    }

    /// Equilibrate a sealed rigid chamber at an explicit total internal energy.
    ///
    /// Use original species formation references when summing feed energies.
    /// No reaction heat is added separately. Source/trace/domain failures and
    /// cancellation return no partially resolved parcel.
    #[allow(clippy::too_many_arguments)]
    pub fn equilibrium_chamber_uv_energy(
        &self,
        chemistry: &IdealGasEquilibriumV1,
        context: &QueryPoint,
        feed: &[(SpeciesId, Amount)],
        target: Energy,
        volume: Volume,
        controls: GasIsochoricEquilibriumSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasIsochoricEquilibriumChamber, GasIsochoricEquilibriumError> {
        self.chamber_uv(
            chemistry, context, feed, target, volume, None, controls, cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn chamber_uv(
        &self,
        chemistry: &IdealGasEquilibriumV1,
        context: &QueryPoint,
        feed: &[(SpeciesId, Amount)],
        target: Energy,
        volume: Volume,
        conserved_mass: Option<f64>,
        controls: GasIsochoricEquilibriumSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasIsochoricEquilibriumChamber, GasIsochoricEquilibriumError> {
        self.admit_equilibrium_bank(chemistry, controls.mass_tolerance_kg, &mut cancelled)?;
        let equilibrium = chemistry.equilibrate_uv_energy(
            volume,
            feed,
            target,
            controls.equilibrium,
            &mut cancelled,
        )?;
        let chemical = equilibrium.equilibrium();
        let feed_mass = equilibrium_feed_mass(chemistry, feed)?;
        let mass = conserved_mass.unwrap_or(feed_mass);
        if !positive(mass) || !positive(mass / volume.value()) {
            return Err(GasIsochoricEquilibriumError::Invalid(
                "unrepresentable conserved density",
            ));
        }
        let (model, correction) = self.freeze_equilibrium_composition(chemical)?;
        // Close the common EOS at the EXACT conserved mass/volume. The chemistry
        // root has a declared residual; that must not accumulate as mass drift.
        let properties = model
            .thermodynamics()
            .evaluate_properties(chemical.temperature(), chemical.pressure())
            .map_err(|error| GasEquilibriumChamberError::from(failure("UV products", error)))?;
        let pressure = (mass / volume.value())
            * properties.mass_specific_gas_constant().value()
            * chemical.temperature().value();
        let point = updated_point(context, chemical.temperature().value(), pressure)
            .map_err(GasEquilibriumChamberError::from)?;
        poll(&mut cancelled)?;
        let state = model
            .evaluate(&point)
            .map_err(GasEquilibriumChamberError::from)?;
        let internal_energy = mass * state.specific_internal_energy_j_kg();
        let balance = GasIsochoricEquilibriumBalance {
            feed_mass_kg: feed_mass,
            mass_residual_kg: mass - feed_mass,
            chemical_mass_residual_kg: chemical.mass().value() - mass,
            density_mass_residual_kg: state.state().density * volume.value() - mass,
            parcel_internal_energy_j: internal_energy,
            energy_residual_j: internal_energy - target.value(),
            freeze_energy_residual_j: internal_energy - equilibrium.internal_energy().value(),
            relative_pressure_correction: pressure / chemical.pressure().value() - 1.0,
            maximum_mole_fraction_correction: correction,
        };
        let pressure_tolerance = 2.0
            * (controls.equilibrium.volume.relative_volume_tolerance
                + controls.mass_tolerance_kg / mass);
        if !internal_energy.is_finite()
            || !balance.energy_residual_j.is_finite()
            || !balance.freeze_energy_residual_j.is_finite()
            || !balance.relative_pressure_correction.is_finite()
            || balance.mass_residual_kg.abs() > controls.mass_tolerance_kg
            || balance.chemical_mass_residual_kg.abs() > controls.mass_tolerance_kg
            || balance.density_mass_residual_kg.abs() > controls.mass_tolerance_kg
            || balance.energy_residual_j.abs() > controls.equilibrium.energy_tolerance_j
            || balance.freeze_energy_residual_j.abs() > controls.equilibrium.energy_tolerance_j
            || balance.relative_pressure_correction.abs() > pressure_tolerance
        {
            return Err(GasIsochoricEquilibriumError::Balance(balance));
        }
        let parcel = ClosedGasParcel::from_isochoric(IsochoricGasParcel {
            model,
            point,
            state,
            mass_kg: mass,
            volume_m3: volume.value(),
        });
        poll(&mut cancelled)?;
        Ok(GasIsochoricEquilibriumChamber {
            equilibrium,
            parcel,
            balance,
        })
    }
}

impl ClosedGasParcel {
    /// Relax a sealed rigid parcel to equilibrium at its actual stored U and V.
    ///
    /// The complete product bank and chemistry must retain the exact active
    /// feed laws and transport sources. Species absent initially may be created
    /// only through the admitted elemental constraints. Mass and volume retain
    /// their exact original bits. The resulting material composition is frozen;
    /// a later kinetic or equilibrium evolution requires another explicit call.
    /// No time increment or instantaneous ignition claim is inferred.
    pub fn equilibrate_adiabatic(
        &mut self,
        bank: &ThermochemicalGasModel,
        chemistry: &IdealGasEquilibriumV1,
        controls: GasIsochoricEquilibriumSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        (
            IdealGasIsochoricEquilibriumV1,
            GasIsochoricEquilibriumBalance,
        ),
        GasIsochoricEquilibriumError,
    > {
        poll(&mut cancelled)?;
        let active = self.model().thermodynamics().components();
        let bank_species = bank.thermodynamics().components();
        for (i, species) in active.iter().enumerate() {
            let Some(j) = bank_species
                .iter()
                .position(|candidate| candidate.species() == species.species())
            else {
                return Err(
                    GasEquilibriumChamberError::SourceMismatch(species.species().clone()).into(),
                );
            };
            if species != &bank_species[j]
                || self.model().transport_materials()[i] != bank.transport_materials()[j]
            {
                return Err(
                    GasEquilibriumChamberError::SourceMismatch(species.species().clone()).into(),
                );
            }
        }
        // The explicit water identity affects vapor fractions and later phase
        // admission even when every numeric source card is identical. A water
        // product absent from the feed may enter through the full bank.
        let active_water = bank
            .water_species()
            .filter(|water| active.iter().any(|species| species.species() == *water));
        if active_water != self.model().water_species() {
            return Err(GasIsochoricEquilibriumError::Invalid(
                "feed and product bank active water identities differ",
            ));
        }
        let active_packs: Vec<_> = bank
            .source_packs()
            .iter()
            .filter(|receipt| {
                active
                    .iter()
                    .any(|species| species.species() == &receipt.species)
            })
            .cloned()
            .collect();
        if active_packs != self.model().source_packs() {
            return Err(GasIsochoricEquilibriumError::Invalid(
                "feed and product bank pack identities differ",
            ));
        }
        let amount = self.mass_kg()
            / self
                .state()
                .thermodynamics()
                .properties()
                .molar_mass()
                .value();
        let feed: Vec<_> = active
            .iter()
            .zip(self.model().thermodynamics().mole_composition().fractions())
            .map(|(species, &fraction)| (species.species().clone(), Amount::new(amount * fraction)))
            .collect();
        if feed.iter().any(|(_, n)| !positive(n.value())) {
            return Err(GasIsochoricEquilibriumError::Invalid(
                "unrepresentable positive feed amount",
            ));
        }
        let target = self.mass_kg() * self.state().specific_internal_energy_j_kg();
        let candidate = bank.chamber_uv(
            chemistry,
            self.point(),
            &feed,
            Energy::new(target),
            Volume::new(self.volume_m3()),
            Some(self.mass_kg()),
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        let (parcel, equilibrium, balance) = candidate.into_parts();
        *self = parcel;
        Ok((equilibrium, balance))
    }
}

#[cfg(test)]
mod tests;
