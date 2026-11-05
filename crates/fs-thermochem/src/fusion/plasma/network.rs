//! Original thermal sources for the bounded light-ion fusion reaction network.

use super::{
    Amount, DtFusionPlasmaErrorV1, DtFusionPlasmaLimitsV1, DtFusionPlasmaScalesV1,
    FrozenIonizedGasStateV1, IdealGasEquilibriumV1, NuclearSpeciesV1, RELATIVE_CHARGE_TOLERANCE,
    SpeciesId, Temperature, Volume, admit_bound_state, admit_inventory, admit_physical_point,
    admit_source_interval, bind_sources,
};

/// Immutable p+/D+/T+/He3++/He4++/electron source binding for light-ion fusion.
///
/// The six original species supply monatomic sensible energy, heat capacities
/// and pressure. Isotopes are explicitly associated with bare-nucleus masses;
/// their names carry no isotope authority. Nuclear rest energy and escaped
/// neutrons remain outside this thermal bank. All heavy populations share one
/// Maxwellian temperature, separate from the electron temperature.
///
/// This binding admits sources and ideal-plasma domains without selecting
/// channels or evolving a burn. The caller selects the original reaction laws
/// and product thermalization assumptions. No chemical element-conservation
/// or equilibrium claim is made for nuclear transmutation.
#[derive(Clone, Debug, PartialEq)]
pub struct LightIonFusionPlasmaV1 {
    thermodynamics: IdealGasEquilibriumV1,
    indices: [usize; 5],
    electron_index: usize,
    log_masses_kg: [f64; 6],
}

impl LightIonFusionPlasmaV1 {
    /// Bind exactly five distinct charged nuclei and one structural electron.
    ///
    /// Original source masses must match the mapped nominal bare nuclei and
    /// electron within 64 epsilon relative tolerance. Every source requires
    /// fixed Cp/R=5/2, zero formation constant a7, fixed a8 across regions and
    /// an untransformed caloric reference. Sources and provenance are retained
    /// unchanged; no replacement ideal-gas cards are constructed.
    pub fn new(
        thermodynamics: IdealGasEquilibriumV1,
        proton: SpeciesId,
        deuterium: SpeciesId,
        tritium: SpeciesId,
        helium3: SpeciesId,
        helium4: SpeciesId,
    ) -> Result<Self, DtFusionPlasmaErrorV1> {
        let (indices, electron_index, log_masses_kg) = bind_sources::<5, 6>(
            &thermodynamics,
            [proton, deuterium, tritium, helium3, helium4],
            [
                NuclearSpeciesV1::Proton,
                NuclearSpeciesV1::Deuterium,
                NuclearSpeciesV1::Tritium,
                NuclearSpeciesV1::Helium3,
                NuclearSpeciesV1::Helium4,
            ],
            "require exactly p+, D+, T+, He3++, He4++, electron on H/He axes",
            "p, D, T, He3 and He4 require distinct source species",
            "fusion mapping requires three H+, two He++ and a structural electron",
        )?;
        Ok(Self {
            thermodynamics,
            indices,
            electron_index,
            log_masses_kg,
        })
    }

    /// Original common source bank, including unchanged source identities.
    pub const fn thermodynamics(&self) -> &IdealGasEquilibriumV1 {
        &self.thermodynamics
    }

    /// Canonical bank indices in the physical order [p, D, T, He3, He4].
    pub const fn indices(&self) -> [usize; 5] {
        self.indices
    }

    /// Canonical bank index of the original nucleus-free charge-minus-one row.
    pub const fn electron_index(&self) -> usize {
        self.electron_index
    }

    /// Check exact source-bank identity and every possible product's support,
    /// including sources whose actual population is still zero.
    pub fn admit_state(
        &self,
        state: &FrozenIonizedGasStateV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), DtFusionPlasmaErrorV1> {
        admit_bound_state(
            &self.thermodynamics,
            self.electron_index,
            state,
            &mut cancelled,
        )
    }

    /// Admit the original thermal state with explicit ideal-plasma limits.
    pub fn admit_state_with_limits(
        &self,
        state: &FrozenIonizedGasStateV1,
        limits: DtFusionPlasmaLimitsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DtFusionPlasmaScalesV1, DtFusionPlasmaErrorV1> {
        self.admit_state(state, &mut cancelled)?;
        self.admit_point(
            state.heavy_temperature(),
            state.electron_temperature(),
            state.amounts(),
            state.volume(),
            limits,
            cancelled,
        )
    }

    /// Allocation-free source, charge and physical-domain admission.
    ///
    /// Amounts follow the retained bank's canonical ordering. Both helium
    /// isotopes carry charge +2 and contribute charge-square four to coupling.
    /// The electron inventory must agree with actual positive carrier charge
    /// to 64 epsilon relative tolerance. No seed, renormalization or population
    /// correction is inserted, including when fuel or secondary products vanish.
    pub fn admit_point(
        &self,
        heavy_temperature: Temperature,
        electron_temperature: Temperature,
        amounts: &[Amount],
        volume: Volume,
        limits: DtFusionPlasmaLimitsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DtFusionPlasmaScalesV1, DtFusionPlasmaErrorV1> {
        let scale = admit_inventory::<6>(
            amounts,
            volume,
            limits,
            "six nonnegative finite amounts and positive finite volume required",
            &mut cancelled,
        )?;
        let [p, d, t, he3, he4] = self.indices;
        let positive_charge = amounts[p].value() / scale
            + amounts[d].value() / scale
            + amounts[t].value() / scale
            + 2.0 * (amounts[he3].value() / scale)
            + 2.0 * (amounts[he4].value() / scale);
        let electron_charge = amounts[self.electron_index].value() / scale;
        if ((positive_charge - electron_charge) / (positive_charge + electron_charge)).abs()
            > RELATIVE_CHARGE_TOLERANCE
        {
            return Err(DtFusionPlasmaErrorV1::Invalid(
                "light-ion trial charge must match the original quasineutral electron inventory",
            ));
        }
        admit_physical_point(
            &self.thermodynamics,
            self.electron_index,
            &self.log_masses_kg,
            heavy_temperature,
            electron_temperature,
            amounts,
            volume,
            limits,
            &mut cancelled,
        )
    }

    /// Admit contiguous original support for each thermal pool, including all
    /// initially absent secondary reactants and products. Original gaps and
    /// unsupported temperatures refuse without extrapolation.
    pub fn admit_temperature_interval(
        &self,
        minimum_heavy: Temperature,
        maximum_heavy: Temperature,
        minimum_electron: Temperature,
        maximum_electron: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), DtFusionPlasmaErrorV1> {
        admit_source_interval(
            &self.thermodynamics,
            self.electron_index,
            minimum_heavy,
            maximum_heavy,
            minimum_electron,
            maximum_electron,
            &mut cancelled,
        )
    }
}
