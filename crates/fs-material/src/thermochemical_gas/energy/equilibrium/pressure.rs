//! Source-bound reacting pressure/enthalpy endpoints for material consumers.
//!
//! The native prepared chemistry path owns composition and its reacting
//! derivatives. The existing material bank freezes only the accepted products,
//! preserving source cards, transport choices and the original feed mass.
//! Material enthalpy increments integrate the original sources directly before
//! rounding; large formation datums do not redefine a small requested change.
//! This supplies thermodynamic endpoints, not reaction kinetics, a throttle
//! rate, a pressure-controller trajectory or an unaccounted heat/work law.

/// Fixed-pressure reacting calorics on the existing finite thermal-storage owner.
pub mod storage;

use super::{GasEquilibriumChamberError, ThermochemicalGasModel};
use crate::thermochemical_gas::{coordinates, energy::parcel::ClosedGasParcel, positive};
use fs_matdb::QueryPoint;
use fs_math::dd::Dd;
use fs_qty::{Amount, Energy, Pressure, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::composition::{
    CompositionEquilibriumErrorV1, EquilibriumEnthalpySolveV1, EquilibriumExpansionErrorV1,
    EquilibriumThermodynamicStateV1, IdealGasEquilibriumV1,
    enthalpy::{EquilibriumEnthalpyPressurePathV1, EquilibriumEnthalpyPressureSampleV1},
};

/// Native work bounds and independent extensive material-balance allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasEquilibriumEnthalpySolve {
    /// Original source-temperature bracket, reacting HP inverse and native H tolerance.
    pub equilibrium: EquilibriumEnthalpySolveV1,
    /// Allowed chemical-minus-feed and density-volume-minus-feed mass [kg].
    pub mass_tolerance_kg: f64,
    /// Allowed material enthalpy residual and endpoint freezing defect [J].
    /// This is an absolute allowance, never scaled by a formation-energy datum.
    pub enthalpy_tolerance_j: f64,
}

/// Actual frozen-material balances against the immutable original reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasEquilibriumEnthalpyBalance {
    /// Exactly preserved original feed mass [kg].
    pub mass_kg: f64,
    /// Native chemical product mass minus the preserved feed mass [kg].
    pub chemical_mass_residual_kg: f64,
    /// Resolved density times accepted volume minus preserved mass [kg].
    pub density_mass_residual_kg: f64,
    /// Actual material H-H_reference from original source integrals [J].
    pub enthalpy_change_j: f64,
    /// Actual material enthalpy change minus the requested change [J].
    pub enthalpy_residual_j: f64,
    /// Direct source H difference between frozen and published chemical amounts [J].
    pub freeze_enthalpy_residual_j: f64,
    /// Rounded absolute material enthalpy [J]; not used to recover increments.
    pub parcel_enthalpy_j: f64,
    /// Largest mole-fraction correction from the existing exact-unit-sum closure.
    pub maximum_mole_fraction_correction: f64,
}

/// Accepted reacting state alongside its original-owner frozen material parcel.
#[derive(Clone, Debug, PartialEq)]
pub struct GasEquilibriumEnthalpyState {
    equilibrium: EquilibriumThermodynamicStateV1,
    parcel: ClosedGasParcel,
    sample: EquilibriumEnthalpyPressureSampleV1,
    balance: GasEquilibriumEnthalpyBalance,
}
impl GasEquilibriumEnthalpyState {
    /// Original chemical source receipts, entropy and reacting response.
    pub const fn equilibrium(&self) -> &EquilibriumThermodynamicStateV1 {
        &self.equilibrium
    }
    /// Frozen products for existing storage, transport, acoustic and flow APIs.
    /// Their frozen Cp and sound speed are distinct from the reacting response.
    pub const fn parcel(&self) -> &ClosedGasParcel {
        &self.parcel
    }
    /// Original native PH residual, reacting derivatives and chemical entropy change.
    pub const fn native_sample(&self) -> EquilibriumEnthalpyPressureSampleV1 {
        self.sample
    }
    /// Actual material mass, enthalpy and composition-freezing discrepancies.
    pub const fn balance(&self) -> GasEquilibriumEnthalpyBalance {
        self.balance
    }
    /// Transfer products without recomputing chemistry or material sources.
    pub fn into_parts(
        self,
    ) -> (
        ClosedGasParcel,
        EquilibriumThermodynamicStateV1,
        EquilibriumEnthalpyPressureSampleV1,
        GasEquilibriumEnthalpyBalance,
    ) {
        (self.parcel, self.equilibrium, self.sample, self.balance)
    }
}

/// No failed query publishes an incomplete or stale material endpoint.
#[derive(Clone, Debug, PartialEq)]
pub enum GasEquilibriumEnthalpyError {
    /// Original source-path, chemistry, reacting derivative or HP refusal.
    Chemistry(EquilibriumExpansionErrorV1),
    /// Original material, neutral transport, source identity or trace refusal.
    Material(GasEquilibriumChamberError),
    /// Invalid or unrepresentable physical data.
    Invalid(&'static str),
    /// Actual material mass, freezing or source-relative H failed its allowance.
    Balance(GasEquilibriumEnthalpyBalance),
    /// Cancellation before publication.
    Cancelled,
}
impl core::fmt::Display for GasEquilibriumEnthalpyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting material pressure/enthalpy path: {self:?}")
    }
}
impl std::error::Error for GasEquilibriumEnthalpyError {}
impl From<EquilibriumExpansionErrorV1> for GasEquilibriumEnthalpyError {
    fn from(error: EquilibriumExpansionErrorV1) -> Self {
        match error {
            EquilibriumExpansionErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled) => {
                Self::Cancelled
            }
            other => Self::Chemistry(other),
        }
    }
}
impl From<GasEquilibriumChamberError> for GasEquilibriumEnthalpyError {
    fn from(error: GasEquilibriumChamberError) -> Self {
        match error {
            GasEquilibriumChamberError::Cancelled => Self::Cancelled,
            other => Self::Material(other),
        }
    }
}
type Error = GasEquilibriumEnthalpyError;
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

/// Reusable chemistry with one immutable material reference and exact feed mass.
///
/// The complete bank retains every potentially formed species. A query updates
/// only numerical warm starts until chemistry, original transport and actual
/// material balances have all passed. Failures invalidate `current`; the
/// reference and all previously returned owned states remain unchanged.
pub struct EquilibriumGasEnthalpyPath<'a> {
    bank: &'a ThermochemicalGasModel,
    context: QueryPoint,
    path: EquilibriumEnthalpyPressurePathV1<'a>,
    controls: GasEquilibriumEnthalpySolve,
    reference: GasEquilibriumEnthalpyState,
    current: Option<GasEquilibriumEnthalpyState>,
}
impl<'a> EquilibriumGasEnthalpyPath<'a> {
    /// Equilibrate the feed elements at the context's original typed T,p.
    /// This initial equilibrium, not unreacted feed enthalpy, defines H_reference.
    pub fn new(
        bank: &'a ThermochemicalGasModel,
        chemistry: &'a IdealGasEquilibriumV1,
        context: &QueryPoint,
        feed: &[(SpeciesId, Amount)],
        controls: GasEquilibriumEnthalpySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if !positive(controls.enthalpy_tolerance_j) {
            return Err(Error::Invalid(
                "positive finite material enthalpy allowance required",
            ));
        }
        bank.admit_equilibrium_bank(chemistry, controls.mass_tolerance_kg, &mut cancelled)?;
        let (temperature, pressure) =
            coordinates(context).map_err(GasEquilibriumChamberError::from)?;
        let mut path = EquilibriumEnthalpyPressurePathV1::new(
            chemistry,
            temperature,
            pressure,
            feed,
            controls.equilibrium,
            &mut cancelled,
        )?;
        let sample = path.sample(pressure, Energy::new(0.0), &mut cancelled)?;
        let native = path.reference_state().clone();
        let reference = freeze(
            bank,
            context,
            native,
            sample,
            path.feed_mass_kg(),
            None,
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            bank,
            context: context.clone(),
            path,
            controls,
            current: Some(reference.clone()),
            reference,
        })
    }

    /// Immutable original material endpoint and chemical source reference.
    pub const fn reference_state(&self) -> &GasEquilibriumEnthalpyState {
        &self.reference
    }

    /// Last complete material query; `None` after any failed or cancelled query.
    pub const fn current(&self) -> Option<&GasEquilibriumEnthalpyState> {
        self.current.as_ref()
    }

    /// Solve p and H-H_reference, then resolve the original material products.
    ///
    /// At constant pressure with only quasistatic pV work the requested change
    /// can equal supplied heat. For other processes the caller must close its
    /// own heat, work and kinetic-energy balance. Native chemical entropy is
    /// reported separately; no reservoir entropy or passivity claim is inferred.
    pub fn sample(
        &mut self,
        pressure: Pressure,
        enthalpy_change: Energy,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&GasEquilibriumEnthalpyState, Error> {
        self.current = None;
        poll(&mut cancelled)?;
        let sample = self
            .path
            .sample(pressure, enthalpy_change, &mut cancelled)?;
        let native = self.path.publish_current(&mut cancelled)?;
        let state = freeze(
            self.bank,
            &self.context,
            native,
            sample,
            self.path.feed_mass_kg(),
            Some(&self.reference.parcel),
            self.controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(self.current.insert(state))
    }
}

#[allow(clippy::too_many_arguments)]
fn freeze(
    bank: &ThermochemicalGasModel,
    context: &QueryPoint,
    equilibrium: EquilibriumThermodynamicStateV1,
    sample: EquilibriumEnthalpyPressureSampleV1,
    mass: f64,
    reference: Option<&ClosedGasParcel>,
    controls: GasEquilibriumEnthalpySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GasEquilibriumEnthalpyState, Error> {
    let chemical = equilibrium.equilibrium();
    let (parcel, correction) =
        bank.freeze_pressure_products(context, chemical, Some(mass), cancelled)?;
    let mixture = parcel.model().thermodynamics();
    let total = Dd::from_f64(mass)
        / Dd::from_f64(
            parcel
                .state()
                .thermodynamics()
                .properties()
                .molar_mass()
                .value(),
        );
    let mut defect = Dd::ZERO;
    for (model, &fraction) in mixture
        .components()
        .iter()
        .zip(mixture.mole_composition().fractions())
    {
        poll(cancelled)?;
        let product = chemical
            .species()
            .iter()
            .find(|s| &s.species == model.species())
            .ok_or(Error::Invalid(
                "frozen product lacks its original chemical species",
            ))?;
        let source = model
            .evaluate_properties(chemical.temperature())
            .map_err(|source| {
                EquilibriumExpansionErrorV1::from(
                    CompositionEquilibriumErrorV1::SpeciesEvaluation {
                        species: model.species().clone(),
                        source,
                    },
                )
            })?;
        defect = defect
            + (total * Dd::from_f64(fraction) - Dd::from_f64(product.amount.value()))
                * Dd::from_f64(source.h().value());
    }
    let enthalpy_change_j = match reference {
        Some(original) => material_enthalpy_change(original, &parcel, cancelled)?,
        None => 0.0,
    };
    let balance = GasEquilibriumEnthalpyBalance {
        mass_kg: mass,
        chemical_mass_residual_kg: chemical.mass().value() - mass,
        density_mass_residual_kg: parcel.state().state().density * parcel.volume_m3() - mass,
        enthalpy_change_j,
        enthalpy_residual_j: enthalpy_change_j - sample.target_enthalpy_change_j,
        freeze_enthalpy_residual_j: defect.to_f64(),
        parcel_enthalpy_j: mass * parcel.state().specific_enthalpy_j_kg(),
        maximum_mole_fraction_correction: correction,
    };
    if [
        balance.chemical_mass_residual_kg,
        balance.density_mass_residual_kg,
        balance.enthalpy_change_j,
        balance.enthalpy_residual_j,
        balance.freeze_enthalpy_residual_j,
        balance.parcel_enthalpy_j,
    ]
    .iter()
    .any(|v| !v.is_finite())
        || balance.chemical_mass_residual_kg.abs() > controls.mass_tolerance_kg
        || balance.density_mass_residual_kg.abs() > controls.mass_tolerance_kg
        || balance.enthalpy_residual_j.abs() > controls.enthalpy_tolerance_j
        || balance.freeze_enthalpy_residual_j.abs() > controls.enthalpy_tolerance_j
    {
        return Err(Error::Balance(balance));
    }
    poll(cancelled)?;
    Ok(GasEquilibriumEnthalpyState {
        equilibrium,
        parcel,
        sample,
        balance,
    })
}

// Independently audit the actual frozen material amounts, including their
// roundoff closure. Published chemical amounts have already been rounded, so
// differencing their absolute total H would not close this material increment.
fn material_enthalpy_change(
    reference: &ClosedGasParcel,
    current: &ClosedGasParcel,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let initial = reference.model().thermodynamics();
    let final_model = current.model().thermodynamics();
    if initial.components() != final_model.components()
        || reference.mass_kg().to_bits() != current.mass_kg().to_bits()
    {
        return Err(Error::Invalid(
            "material PH reference must retain original active sources and mass",
        ));
    }
    let initial_total = Dd::from_f64(reference.mass_kg())
        / Dd::from_f64(
            reference
                .state()
                .thermodynamics()
                .properties()
                .molar_mass()
                .value(),
        );
    let final_total = Dd::from_f64(current.mass_kg())
        / Dd::from_f64(
            current
                .state()
                .thermodynamics()
                .properties()
                .molar_mass()
                .value(),
        );
    let initial_t = Temperature::new(reference.state().state().temperature);
    let final_t = Temperature::new(current.state().state().temperature);
    let mut change = Dd::ZERO;
    for ((model, &before), &after) in initial
        .components()
        .iter()
        .zip(initial.mole_composition().fractions())
        .zip(final_model.mole_composition().fractions())
    {
        poll(cancelled)?;
        let source_error = |source| {
            EquilibriumExpansionErrorV1::from(CompositionEquilibriumErrorV1::SpeciesEvaluation {
                species: model.species().clone(),
                source,
            })
        };
        let h = model
            .evaluate_properties(initial_t)
            .map_err(source_error)?
            .h()
            .value();
        let sensible = model
            .molar_enthalpy_change(initial_t, final_t)
            .map_err(source_error)?
            .value();
        let n = final_total * Dd::from_f64(after);
        let dn = n - initial_total * Dd::from_f64(before);
        change = change + dn * Dd::from_f64(h) + n * Dd::from_f64(sensible);
    }
    let change = change.to_f64();
    if !change.is_finite() {
        return Err(Error::Invalid(
            "material enthalpy increment is unrepresentable",
        ));
    }
    Ok(change)
}

#[cfg(test)]
mod tests;
