//! Sealed reacting gas on the original finite thermal-storage owner.
//!
//! Cells store cold-relative internal energy at fixed mass/volume ratio and
//! elemental composition. The original TV/UV closure supplies reacting Cv and
//! evolving pressure; the existing thermal timeline owns every transferred
//! joule. There is no additional reaction-heat source or pressure-volume work
//! for this rigid constraint. This is the infinitely fast neutral ideal-gas
//! equilibrium limit, not a finite-rate ignition or open-mass flow model.

/// Original reacting calorics with an explicitly declared bulk k(T) curve.
pub mod conduction;

use super::{
    GasEquilibriumChamberError, GasIsochoricEquilibriumChamber, GasIsochoricEquilibriumError,
    ThermochemicalGasModel,
    driven::{DrivenEquilibriumGasError as Error, DrivenEquilibriumGasSolve, freeze},
};
use crate::thermochemical_gas::{energy::GasEnergySolve, positive};
use fs_matdb::QueryPoint;
use fs_math::dd::Dd;
use fs_qty::{Amount, Temperature, Volume};
use fs_thermal::enthalpy::{DepositedCell, DepositionError, ThermalLaw, ThermalStorageBasis};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::composition::{
    IdealGasEquilibriumV1,
    volume::energy_path::{
        EquilibriumEnergyVolumeSampleV1, isochoric::EquilibriumIsochoricCaloricV1,
    },
};
use std::cell::RefCell;

/// Accepted original material endpoint and its source-relative storage balance.
#[derive(Clone, Debug, PartialEq)]
pub struct GasIsochoricThermalState {
    chamber: GasIsochoricEquilibriumChamber,
    sample: EquilibriumEnergyVolumeSampleV1,
    internal_energy_change_j: f64,
    energy_residual_j: f64,
    entropy_change_j_k: f64,
}
impl GasIsochoricThermalState {
    /// Original native chemistry and separately frozen material transport state.
    pub const fn chamber(&self) -> &GasIsochoricEquilibriumChamber {
        &self.chamber
    }
    /// Original UV response, energy/volume residuals and bounded work counts.
    pub const fn native_sample(&self) -> EquilibriumEnergyVolumeSampleV1 {
        self.sample
    }
    /// Actual material U-U_cold, including original freezing corrections [J].
    pub const fn internal_energy_change_j(&self) -> f64 {
        self.internal_energy_change_j
    }
    /// Actual material change minus the requested stored-energy change [J].
    pub const fn energy_residual_j(&self) -> f64 {
        self.energy_residual_j
    }
    /// Actual material S-S_cold, including original freezing corrections [J/K].
    pub const fn entropy_change_j_k(&self) -> f64 {
        self.entropy_change_j_k
    }
}

/// Immutable physical calorics with bounded, canonically restored solver scratch.
///
/// A thermal cell may use any positive mass of the prepared element ratios,
/// provided its rigid volume retains the prepared density. Open-mass consumers
/// see an explicit InternalEnergy basis and no inferred advection descriptor.
pub struct GasIsochoricEquilibriumThermalStorage<'a> {
    bank: &'a ThermochemicalGasModel,
    chemistry: &'a IdealGasEquilibriumV1,
    context: QueryPoint,
    feed: Vec<(SpeciesId, Amount)>,
    controls: DrivenEquilibriumGasSolve,
    reference: GasIsochoricEquilibriumChamber,
    reference_freeze_energy_j: f64,
    reference_freeze_entropy_j_k: f64,
    caloric: RefCell<EquilibriumIsochoricCaloricV1<'a>>,
    mass_kg: f64,
    density_kg_m3: f64,
    temperature_domain: [Temperature; 2],
    energy_domain_j_kg: [f64; 2],
    maximum_slope: f64,
}
impl<'a> GasIsochoricEquilibriumThermalStorage<'a> {
    /// Prepare the original continuous source path and material endpoints.
    ///
    /// The volume and feed, not the context pressure, fix density. Auxiliary
    /// material axes are retained; original product freezing re-closes pressure
    /// at exact mass/volume and checks its admitted correction. Caloric trials
    /// use complete NASA support; accepted full material queries separately
    /// recheck transport claims at the evolving pressure and composition.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bank: &'a ThermochemicalGasModel,
        chemistry: &'a IdealGasEquilibriumV1,
        context: &QueryPoint,
        feed: &[(SpeciesId, Amount)],
        volume: Volume,
        controls: DrivenEquilibriumGasSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if !positive(controls.energy_tolerance_j) || !positive(controls.entropy_tolerance_j_k) {
            return Err(Error::Invalid(
                "positive finite material U and entropy allowances required",
            ));
        }
        bank.admit_equilibrium_bank(chemistry, controls.mass_tolerance_kg, &mut cancelled)
            .map_err(GasIsochoricEquilibriumError::from)?;
        let mut caloric = EquilibriumIsochoricCaloricV1::new(
            chemistry,
            volume,
            feed,
            controls.equilibrium,
            &mut cancelled,
        )?;
        let mass_kg = caloric.feed_mass_kg();
        let density_kg_m3 = caloric.density_kg_m3();
        let temperature_domain = caloric.temperature_domain();
        let energy_domain_j_kg = caloric.specific_internal_energy_domain_j_kg();
        let cold = caloric.sample_specific_internal_energy(0.0, &mut cancelled)?;
        let (reference, reference_freeze_energy_j, reference_freeze_entropy_j_k) = freeze(
            bank,
            chemistry,
            context,
            feed,
            caloric.publish_current(&mut cancelled)?,
            mass_kg,
            controls,
            &mut cancelled,
        )?;

        // Minimizing ideal-gas Helmholtz energy over the admitted composition
        // adds a nonnegative redistribution term to Cv at fixed V. Therefore
        // Cv_eq >= sum(n_i cv_i). Reuse the original pure-species coefficient
        // bounds over EVERY eligible frozen-reference source, not a frozen
        // mixture capacity sample. A positive native trace that cannot enter
        // the material representation has already refused at the boundary.
        let source_controls = GasEnergySolve {
            minimum_temperature_k: temperature_domain[0].value(),
            maximum_temperature_k: temperature_domain[1].value(),
            energy_tolerance_j_kg: controls.equilibrium.energy_tolerance_j / mass_kg,
            maximum_iterations: controls.equilibrium.maximum_temperature_iterations.min(128),
        };
        let minimum_cv = reference
            .parcel()
            .model()
            .prepare_caloric_branch(source_controls, &mut cancelled)
            .map_err(GasEquilibriumChamberError::from)
            .map_err(GasIsochoricEquilibriumError::from)?
            .minimum_cv_j_kg_k();
        let retained_mass = (1.0
            - controls
                .equilibrium
                .volume
                .equilibrium
                .relative_balance_tolerance)
            .next_down();
        let lower_capacity = (minimum_cv * retained_mass).next_down();
        let maximum_slope = (1.0 / lower_capacity).next_up();
        if !positive(lower_capacity) || !positive(maximum_slope) {
            return Err(Error::Invalid(
                "no represented source Cv bound for reacting storage",
            ));
        }
        bound_sample(cold, mass_kg, maximum_slope)?;
        let hot = caloric.sample_temperature(temperature_domain[1], &mut cancelled)?;
        bound_sample(hot, mass_kg, maximum_slope)?;
        freeze(
            bank,
            chemistry,
            context,
            feed,
            caloric.publish_current(&mut cancelled)?,
            mass_kg,
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        let law = Self {
            bank,
            chemistry,
            context: context.clone(),
            feed: feed.to_vec(),
            controls,
            reference,
            reference_freeze_energy_j,
            reference_freeze_entropy_j_k,
            caloric: RefCell::new(caloric),
            mass_kg,
            density_kg_m3,
            temperature_domain,
            energy_domain_j_kg,
            maximum_slope,
        };
        // Admit the exact advertised inclusive stored-energy endpoint through
        // the same material relative-U gate used by consumers. Separate cold
        // and hot freeze defects may each pass while their difference does not.
        law.resolve_specific_internal_energy(energy_domain_j_kg[1], &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(law)
    }

    /// Original cold material state, exact inventory/volume and source receipts.
    pub const fn reference_state(&self) -> &GasIsochoricEquilibriumChamber {
        &self.reference
    }
    /// Extensive chemistry normalization mass [kg], not a cell-size restriction.
    pub const fn feed_mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Conserved mass divided by the authored rigid volume [kg/m3].
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    /// Complete continuous original-source temperature interval.
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        self.temperature_domain
    }
    /// Inclusive cold-relative stored-u interval [J/kg].
    pub const fn specific_internal_energy_domain_j_kg(&self) -> [f64; 2] {
        self.energy_domain_j_kg
    }

    /// One canonical original UV inverse supplies T, pressure and reacting Cv.
    pub fn caloric_state_at_specific_internal_energy(
        &self,
        u_j_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumEnergyVolumeSampleV1, Error> {
        poll(&mut cancelled)?;
        let mut caloric = self
            .caloric
            .try_borrow_mut()
            .map_err(|_| Error::Invalid("isochoric caloric query is already in progress"))?;
        let sample = caloric.sample_specific_internal_energy(u_j_kg, &mut cancelled)?;
        bound_sample(sample, self.mass_kg, self.maximum_slope)?;
        poll(&mut cancelled)?;
        Ok(sample)
    }

    /// Direct original TV state expressed as (u [J/kg], reacting Cv [J/(kg K)]).
    pub fn internal_energy_and_heat_capacity_at_temperature(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(f64, f64), Error> {
        poll(&mut cancelled)?;
        let mut caloric = self
            .caloric
            .try_borrow_mut()
            .map_err(|_| Error::Invalid("isochoric caloric query is already in progress"))?;
        let sample = caloric.sample_temperature(temperature, &mut cancelled)?;
        bound_sample(sample, self.mass_kg, self.maximum_slope)?;
        let u = sample.internal_energy_change_j / self.mass_kg;
        let cv = sample.response.heat_capacity_v_j_k / self.mass_kg;
        if !self.contains(u) || !positive(cv) {
            return Err(Error::Invalid(
                "unrepresentable reacting isochoric forward state",
            ));
        }
        poll(&mut cancelled)?;
        Ok((u, cv))
    }

    /// Original complete material endpoint at the same stored u. Material
    /// pressure is reclosed at exact mass/volume by the shared freezing owner;
    /// its correction and source/transport/freezing refusals remain visible.
    pub fn resolve_specific_internal_energy(
        &self,
        u_j_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasIsochoricThermalState, Error> {
        poll(&mut cancelled)?;
        let (sample, native) = {
            let mut caloric = self
                .caloric
                .try_borrow_mut()
                .map_err(|_| Error::Invalid("isochoric caloric query is already in progress"))?;
            let sample = caloric.sample_specific_internal_energy(u_j_kg, &mut cancelled)?;
            bound_sample(sample, self.mass_kg, self.maximum_slope)?;
            (sample, caloric.publish_current(&mut cancelled)?)
        };
        let (chamber, freeze_energy, freeze_entropy) = freeze(
            self.bank,
            self.chemistry,
            &self.context,
            &self.feed,
            native,
            self.mass_kg,
            self.controls,
            &mut cancelled,
        )?;
        // Reuse the existing driven parcel's original-source balance, retaining
        // both tails while combining native changes and amount-only freezing.
        let stored = Dd::from_f64(sample.internal_energy_change_j) + Dd::from_f64(freeze_energy)
            - Dd::from_f64(self.reference_freeze_energy_j);
        let internal_energy_change_j = stored.to_f64();
        let energy_residual_j = (stored - Dd::from_f64(sample.target_energy_change_j)).to_f64();
        let entropy_change_j_k = (Dd::from_f64(sample.entropy_change_j_k)
            + Dd::from_f64(freeze_entropy)
            - Dd::from_f64(self.reference_freeze_entropy_j_k))
        .to_f64();
        if !internal_energy_change_j.is_finite()
            || !energy_residual_j.is_finite()
            || !entropy_change_j_k.is_finite()
            || energy_residual_j.abs() > self.controls.energy_tolerance_j
        {
            return Err(Error::Invalid(
                "original material U does not close the stored thermal state",
            ));
        }
        poll(&mut cancelled)?;
        Ok(GasIsochoricThermalState {
            chamber,
            sample,
            internal_energy_change_j,
            energy_residual_j,
            entropy_change_j_k,
        })
    }

    /// Source inverse-temperature mean times actual stored du, with chemical
    /// mixing and the evolving equilibrium pressure included exactly once.
    /// Original UV/volume accuracy limits resolvable intervals. Unclosed local
    /// energy or a mean outside the actual reciprocal-T interval refuses;
    /// neither entropy nor a second-law balance is clamped into acceptance.
    /// Coincident accepted temperatures use the continuous du/T limit, whose
    /// numerical accuracy still depends on the original absolute UV tolerance.
    pub fn thermal_entropy_change_j_kg_k(
        &self,
        initial_u: f64,
        final_u: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        poll(&mut cancelled)?;
        let change = self
            .caloric
            .try_borrow_mut()
            .map_err(|_| Error::Invalid("isochoric caloric query is already in progress"))?
            .change_between_specific_internal_energies(initial_u, final_u, &mut cancelled)?;
        bound_sample(change.initial, self.mass_kg, self.maximum_slope)?;
        bound_sample(change.final_state, self.mass_kg, self.maximum_slope)?;
        if initial_u == final_u {
            poll(&mut cancelled)?;
            return Ok(0.0);
        }
        let du = final_u - initial_u;
        let residual = (Dd::from_f64(change.internal_energy_change_j)
            - Dd::from_f64(change.requested_internal_energy_change_j)
            - Dd::from_f64(change.final_state.energy_residual_j)
            + Dd::from_f64(change.initial.energy_residual_j))
        .to_f64();
        if !residual.is_finite() || residual.abs() > self.controls.equilibrium.energy_tolerance_j {
            return Err(Error::Invalid(
                "reacting local source U does not close its thermal interval",
            ));
        }
        let (t0, t1) = (
            change.initial.temperature.value(),
            change.final_state.temperature.value(),
        );
        if (du > 0.0 && t1 < t0) || (du < 0.0 && t1 > t0) {
            return Err(Error::Invalid(
                "reacting isochoric inverse is not monotone on the interval",
            ));
        }
        let mean = if t0 == t1 {
            1.0 / t0
        } else {
            if change.internal_energy_change_j == 0.0
                || change.entropy_change_j_k == 0.0
                || change.internal_energy_change_j.is_sign_positive() != du.is_sign_positive()
                || change.entropy_change_j_k.is_sign_positive() != du.is_sign_positive()
            {
                return Err(Error::Invalid(
                    "reacting isochoric source interval is unresolved",
                ));
            }
            let mean = change.entropy_change_j_k / change.internal_energy_change_j;
            let lower = (1.0 / t0.max(t1)).next_down();
            let upper = (1.0 / t0.min(t1)).next_up();
            let allowance = (64.0 * f64::EPSILON * upper).next_up();
            let violation = if mean < lower {
                lower - mean
            } else if mean > upper {
                mean - upper
            } else {
                0.0
            };
            if !positive(mean)
                || !positive(lower)
                || !positive(upper)
                || !allowance.is_finite()
                || violation > allowance
            {
                return Err(Error::Invalid(
                    "reacting isochoric entropy mean is outside its thermal interval",
                ));
            }
            mean
        };
        let entropy = du * mean;
        if !positive(mean) || !entropy.is_finite() || entropy == 0.0 {
            return Err(Error::Invalid(
                "unrepresentable reacting isochoric entropy increment",
            ));
        }
        poll(&mut cancelled)?;
        Ok(entropy)
    }

    fn contains(&self, u: f64) -> bool {
        u.is_finite() && (self.energy_domain_j_kg[0]..=self.energy_domain_j_kg[1]).contains(&u)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn bound_sample(
    sample: EquilibriumEnergyVolumeSampleV1,
    mass: f64,
    maximum_slope: f64,
) -> Result<(), Error> {
    let slope = mass / sample.response.heat_capacity_v_j_k;
    if !positive(slope) || slope > maximum_slope || !positive(sample.pressure.value()) {
        return Err(Error::Invalid(
            "reacting Cv violates its original source capacity bound",
        ));
    }
    Ok(())
}

impl ThermalLaw for GasIsochoricEquilibriumThermalStorage<'_> {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        let expected = self.density_kg_m3 * cell.geometry.volume_m3;
        if !positive(cell.mass_kg)
            || !positive(cell.geometry.volume_m3)
            || !positive(expected)
            || (cell.mass_kg - expected).abs() > 64.0 * f64::EPSILON * cell.mass_kg.max(expected)
        {
            return Err(DepositionError::InvalidInput(
                "cell mass/volume differs from its reacting rigid gas",
            ));
        }
        self.temperature_k(cell.specific_h_j_kg)?;
        Ok(())
    }
    fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
        self.caloric_state_at_specific_internal_energy(u, || false)
            .map(|sample| sample.temperature.value())
            .map_err(|_| DepositionError::MaterialDomain)
    }
    fn bulk_density_kg_m3(&self, u: f64) -> Result<Option<f64>, DepositionError> {
        // Density is the exact immutable inventory/volume constraint; obtaining
        // it needs domain admission, not a second chemistry/temperature solve.
        if self.contains(u) {
            Ok(Some(self.density_kg_m3))
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_slope
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.thermal_entropy_change_j_kg_k(a, b, || false)
            .map(Some)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

#[cfg(test)]
mod tests;
