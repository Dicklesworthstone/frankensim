//! Sourced thermal expansion around the original permanent-creep owners.
//!
//! The original constructors retain fixed E/nu. Explicit thermal constructors
//! also admit independently affine E(T)/nu(T); energetic hardening stays fixed.
//! Actual temperature controls Arrhenius kinetics and the selected yield offset.
//! The existing affine thermoelastic potential owns integrated expansion,
//! entropy and its full-interval capacity bound. There is no second return map.
//!
//! With e = epsilon - epsilon_c - epsilon_p - f(T) I, the native increment
//! supplies free-energy storage and physical/numerical losses. Reciprocal heat
//! is sigma_1:I delta(f) - delta(Psi_T) - delta(T s), where delta(Psi_T) is the
//! change of elastic storage at the accepted elastic strain due to C(T) alone.
//! Thus total stress work equals change of
//! free energy plus delta(T s), physical heat, reciprocal heat and numerical
//! loss. Only the two physical heat-transfer terms enter the caloric chart.

use std::sync::Arc;

use crate::creep::plastic::NortonJ2State;
use crate::creep::plastic::source::NortonJ2Source;
use crate::creep::recovery::plastic::source::RecoveringNortonJ2Source;
use crate::creep::recovery::source::RecoveringNortonSource;
use crate::creep::source::{
    CreepSourceError, CreepSourceFamily, NortonCreepSource, NortonSourceQuery,
};
use crate::creep::{
    CreepError, CreepIncrement, CreepLaw, CreepState, CreepTolerances, NortonCreep,
};
use crate::elastic::thermal::{
    AffineIsotropicThermoelastic, ResolvedThermoelastic, ThermoelasticError, ThermoelasticState,
};
use crate::plastic::source::HardeningPropertyBasis;
use crate::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, ResolvedMaterialStatePoint, integrate_isotropic_thermal_expansion,
};
use crate::tensor::{contract, trace};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::ModelCard;
use fs_matdb::MaterialCard;

/// Original irreversible strain carried by a native creep history.
/// Every component history remains separately bounded by its owning law.
pub trait CreepHistoryStrain {
    /// Number of independently bounded irreversible strain tensors.
    const STRAIN_FIELDS: usize;
    fn inelastic_strain(&self) -> Voigt;
}
impl CreepHistoryStrain for CreepState {
    const STRAIN_FIELDS: usize = 1;
    fn inelastic_strain(&self) -> Voigt {
        self.creep_strain
    }
}
impl CreepHistoryStrain for NortonJ2State {
    const STRAIN_FIELDS: usize = 2;
    fn inelastic_strain(&self) -> Voigt {
        core::array::from_fn(|i| self.creep.creep_strain[i] + self.plastic.plastic_strain[i])
    }
}

/// Retain the actual accepted temperature alongside every original memory.
/// A thermal trial always borrows this complete history, never ages it in place.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalCreepState<S> {
    pub mechanical: S,
    temperature_k: f64,
}
impl<S> ThermalCreepState<S> {
    /// Reconstruct a retained thermal history without advancing either owner.
    pub fn from_parts(mechanical: S, temperature_k: f64) -> Self {
        Self {
            mechanical,
            temperature_k,
        }
    }
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
}

/// Immutable expansion and mechanical sources shared by all physical trials.
#[derive(Clone, Debug)]
pub struct ExpandingCreepSource<S> {
    mechanical: S,
    potential: Arc<AffineIsotropicThermoelastic>,
}
impl<'a, S: CreepSourceFamily<'a>> ExpandingCreepSource<S>
where
    <S::Material as SmallStrainLaw>::State: CreepHistoryStrain,
{
    /// Add an explicitly selected expansion source and stress-free reference.
    /// This source uses the thermoelastic owner's canonical absolute T axis.
    /// Missing expansion data is an error, never a zero-expansion fallback.
    pub fn new(
        card: &'a MaterialCard,
        mechanical: S,
        expansion_reference: &ResolvedMaterialStatePoint,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CreepSourceError> {
        poll(&mut cancelled)?;
        let source = mechanical.creep_source();
        if source.temperature_axis() != "T" {
            return Err(CreepSourceError::Invalid(
                "the thermoelastic expansion source requires the original canonical T axis",
            ));
        }
        let reference = source.reference();
        let expansion = expansion_reference
            .property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
            .ok_or(CreepSourceError::Invalid(
                "explicit thermal expansion source required",
            ))?;
        if reference.card_identity() != card.content_hash()
            || expansion_reference.card_identity() != card.content_hash()
            || !expansion.has_only_fixed_external_conditions(Some("T"), false)
        {
            return Err(CreepSourceError::Invalid(
                "expansion and creep require one original card and fixed external conditions",
            ));
        }
        // The potential sees epsilon - sum(inelastic strains). Its full-interval
        // capacity envelope includes each original tensor's independent bound;
        // cancellation cannot conceal a larger internal strain.
        let envelope = (1
            + <<S::Material as SmallStrainLaw>::State as CreepHistoryStrain>::STRAIN_FIELDS)
            as f64
            * source.strain_limit();
        let potential = AffineIsotropicThermoelastic::new(
            card,
            reference,
            expansion_reference,
            source.parameters().temperature_range_k,
            envelope,
        )
        .map_err(|e| CreepSourceError::Mechanical(thermal_error(e)))?;
        poll(&mut cancelled)?;
        // This owner admits the complete smooth E/nu/alpha paths and proves that
        // the full thermal-capacity bound, including coefficient work, is finite.
        potential
            .minimum_heat_capacity_j_k(1.0)
            .map_err(|e| CreepSourceError::Mechanical(thermal_error(e)))?;
        Ok(Self {
            mechanical,
            potential: Arc::new(potential),
        })
    }

    pub fn mechanical_source(&self) -> &S {
        &self.mechanical
    }
    pub fn potential(&self) -> &AffineIsotropicThermoelastic {
        &self.potential
    }
    pub fn at_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<ThermalCreep<S::Material>, CreepSourceError> {
        let thermal = self
            .potential
            .at_temperature(temperature_k)
            .map_err(|e| CreepSourceError::Mechanical(thermal_error(e)))?;
        let mechanical = self
            .mechanical
            .at_thermoelastic_temperature(temperature_k, duration_s)?;
        Ok(ThermalCreep {
            mechanical,
            thermal,
            potential: self.potential.clone(),
        })
    }
}

impl<'a> ExpandingCreepSource<NortonCreepSource<'a>> {
    /// Bind sourced Norton kinetics and independently affine E(T), nu(T), alpha(T)
    /// to one thermal potential. Every nonelastic Norton coefficient remains
    /// constant over the admitted span. Missing expansion is never inferred zero.
    pub fn from_norton(
        card: &'a MaterialCard,
        query: NortonSourceQuery<'_>,
        expansion_reference: &ResolvedMaterialStatePoint,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CreepSourceError> {
        let mechanical =
            NortonCreepSource::prepare_thermoelastic(card, query, false, &mut cancelled)?;
        Self::new(card, mechanical, expansion_reference, cancelled)
    }

    /// The same thermal path with three explicitly selected primary-creep
    /// claims. Accumulated irreversible strain remains with the native law.
    pub fn from_norton_with_strain_hardening(
        card: &'a MaterialCard,
        query: NortonSourceQuery<'_>,
        expansion_reference: &ResolvedMaterialStatePoint,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CreepSourceError> {
        let mechanical =
            NortonCreepSource::prepare_thermoelastic(card, query, true, &mut cancelled)?;
        Self::new(card, mechanical, expansion_reference, cancelled)
    }
}

impl<'a> ExpandingCreepSource<NortonJ2Source<'a>> {
    /// Simultaneous sourced creep, yield softening and affine elastic softening.
    /// Energetic hardening must use a constant, directly sourced plastic-strain
    /// modulus; a post-yield tangent would change its converted storage with E.
    pub fn from_norton_j2(
        card: &'a MaterialCard,
        query: NortonSourceQuery<'_>,
        basis: HardeningPropertyBasis,
        expansion_reference: &ResolvedMaterialStatePoint,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CreepSourceError> {
        let mechanical =
            NortonJ2Source::prepare_thermoelastic(card, query, basis, false, &mut cancelled)?;
        Self::new(card, mechanical, expansion_reference, cancelled)
    }

    /// Affine thermoelasticity with independently retained primary-creep and J2
    /// plastic histories. All source roles are explicitly selected and pinned.
    pub fn from_norton_j2_with_strain_hardening(
        card: &'a MaterialCard,
        query: NortonSourceQuery<'_>,
        basis: HardeningPropertyBasis,
        expansion_reference: &ResolvedMaterialStatePoint,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CreepSourceError> {
        let mechanical =
            NortonJ2Source::prepare_thermoelastic(card, query, basis, true, &mut cancelled)?;
        Self::new(card, mechanical, expansion_reference, cancelled)
    }
}

impl<'a> ExpandingCreepSource<RecoveringNortonSource<'a>> {
    /// Affine sourced elasticity/expansion with independently recovering primary
    /// resistance. Flow and recovery retain their separate Arrhenius parameters;
    /// recovery never erases permanent creep or supplies stored energy or heat.
    pub fn from_recovering_norton(
        card: &'a MaterialCard,
        query: NortonSourceQuery<'_>,
        expansion_reference: &ResolvedMaterialStatePoint,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CreepSourceError> {
        let mechanical =
            RecoveringNortonSource::prepare_thermoelastic(card, query, &mut cancelled)?;
        Self::new(card, mechanical, expansion_reference, cancelled)
    }
}

impl<'a> ExpandingCreepSource<RecoveringNortonJ2Source<'a>> {
    /// Simultaneous elastic/yield softening, permanent creep/plasticity and
    /// resistance recovery with one thermal energy owner. Energetic hardening
    /// must remain a constant directly sourced plastic-strain modulus.
    pub fn from_recovering_norton_j2(
        card: &'a MaterialCard,
        query: NortonSourceQuery<'_>,
        basis: HardeningPropertyBasis,
        expansion_reference: &ResolvedMaterialStatePoint,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CreepSourceError> {
        let mechanical =
            RecoveringNortonJ2Source::prepare_thermoelastic(card, query, basis, &mut cancelled)?;
        Self::new(card, mechanical, expansion_reference, cancelled)
    }
}

impl<'a, S: CreepSourceFamily<'a>> CreepSourceFamily<'a> for ExpandingCreepSource<S>
where
    <S::Material as SmallStrainLaw>::State: CreepHistoryStrain,
{
    type Material = ThermalCreep<S::Material>;
    fn creep_source(&self) -> &NortonCreepSource<'a> {
        self.mechanical.creep_source()
    }
    fn at_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<Self::Material, CreepSourceError> {
        ExpandingCreepSource::at_temperature(self, temperature_k, duration_s)
    }
    fn minimum_thermal_capacity_j_m3_k(&self) -> Result<f64, CreepSourceError> {
        self.potential
            .minimum_heat_capacity_j_k(1.0)
            .map_err(|e| CreepSourceError::Mechanical(thermal_error(e)))
    }
    fn admit_stress_free_temperature(&self, temperature_k: f64) -> Result<(), CreepSourceError> {
        let thermal = self
            .potential
            .at_temperature(temperature_k)
            .map_err(|e| CreepSourceError::Mechanical(thermal_error(e)))?;
        if thermal
            .free_strain()
            .map_err(|e| CreepSourceError::Mechanical(thermal_error(e)))?
            != 0.0
        {
            return Err(CreepSourceError::Invalid(
                "virgin zero-displacement creep birth requires its sourced stress-free temperature",
            ));
        }
        Ok(())
    }
    fn validate_accepted_temperature(
        &self,
        temperature_k: f64,
        strain: &Voigt,
        state: &ThermalCreepState<<S::Material as SmallStrainLaw>::State>,
    ) -> Result<(), CreepSourceError> {
        if state.temperature_k() != temperature_k {
            return Err(CreepSourceError::Invalid(
                "accepted creep history and caloric temperature must agree",
            ));
        }
        self.at_temperature(temperature_k, 1.0)?
            .validate_creep_history(strain, state)
            .map_err(CreepSourceError::Mechanical)
    }
}

/// Endpoint-temperature constitutive law. All Newton trials retain the old
/// mechanical state and temperature; accepted state is returned only once.
#[derive(Clone, Debug)]
pub struct ThermalCreep<M> {
    mechanical: M,
    thermal: ResolvedThermoelastic,
    potential: Arc<AffineIsotropicThermoelastic>,
}
impl<M: CreepLaw> ThermalCreep<M>
where
    M::State: CreepHistoryStrain,
{
    pub fn mechanical(&self) -> &M {
        &self.mechanical
    }
    pub fn thermoelastic(&self) -> &ResolvedThermoelastic {
        &self.thermal
    }

    /// Original expansion state retained at this endpoint, including its claim receipt.
    pub fn expansion_material(&self) -> &ResolvedMaterialStatePoint {
        &self.thermal.expansion
    }

    /// Integrated sourced free strain at any admitted history temperature.
    pub fn free_strain_at_temperature(&self, temperature_k: f64) -> Result<f64, CreepError> {
        self.potential
            .at_temperature(temperature_k)
            .map_err(thermal_error)?
            .free_strain()
            .map_err(thermal_error)
    }

    /// Integrated sourced free strain at this law's endpoint temperature.
    pub fn current_free_strain(&self) -> Result<f64, CreepError> {
        self.thermal.free_strain().map_err(thermal_error)
    }

    // Reconstruct an accepted free-sided tensor with its ORIGINAL Poisson ratio.
    // The endpoint law may have a different nu and cannot reconstruct old lateral
    // strain without changing the accepted energy/history behind the caller.
    pub(super) fn uniaxial_history_strain(
        &self,
        axial_strain: f64,
        state: &M::State,
        temperature_k: f64,
    ) -> Result<Voigt, CreepError> {
        let thermal = self
            .potential
            .at_temperature(temperature_k)
            .map_err(thermal_error)?;
        let free = thermal.free_strain().map_err(thermal_error)?;
        let poisson = thermal
            .elastic
            .property(POISSON_RATIO_PROPERTY)
            .ok_or(CreepError::Invalid(
                "original Poisson-ratio claim is missing",
            ))?
            .value_si();
        let inelastic = state.inelastic_strain();
        let transverse_elastic = -poisson * (axial_strain - free - inelastic[0]);
        let strain = [
            axial_strain,
            transverse_elastic + inelastic[1] + free,
            transverse_elastic + inelastic[2] + free,
            inelastic[3],
            inelastic[4],
            inelastic[5],
        ];
        let mechanical = self.mechanical_strain(&strain, &thermal)?;
        self.mechanical.validate_creep_history(&mechanical, state)?;
        self.thermal_state(&strain, state, &thermal)?;
        Ok(strain)
    }

    fn mechanical_strain(
        &self,
        total: &Voigt,
        thermal: &ResolvedThermoelastic,
    ) -> Result<Voigt, CreepError> {
        let limit = self.kinetics().elastic().strain_limit;
        let norm = contract(total, total);
        if total.iter().any(|x| !x.is_finite()) || !norm.is_finite() || norm > limit * limit {
            return Err(CreepError::Invalid(
                "total strain leaves the thermocreep envelope",
            ));
        }
        let free = thermal.free_strain().map_err(thermal_error)?;
        Ok(core::array::from_fn(|i| {
            total[i] - if i < 3 { free } else { 0.0 }
        }))
    }

    fn thermal_state(
        &self,
        total: &Voigt,
        state: &M::State,
        thermal: &ResolvedThermoelastic,
    ) -> Result<ThermoelasticState, CreepError> {
        let inelastic = state.inelastic_strain();
        let mut strain = core::array::from_fn(|i| total[i] - inelastic[i]);
        // Native solid/material shear is tensor Voigt; the potential uses Mandel.
        let root_two = fs_math::det::sqrt(2.0);
        for value in &mut strain[3..] {
            *value *= root_two;
        }
        thermal.state(strain, 1.0).map_err(thermal_error)
    }

    fn trial_strain(
        &self,
        strain: &Voigt,
        state: &ThermalCreepState<M::State>,
    ) -> Result<Voigt, CreepError> {
        let [lower, upper] = self.potential.temperature_range_k();
        if !state.temperature_k.is_finite() || !(lower..=upper).contains(&state.temperature_k) {
            return Err(CreepError::Invalid(
                "accepted thermocreep temperature leaves its source span",
            ));
        }
        self.mechanical_strain(strain, &self.thermal)
    }

    /// Audit the original return map and the reciprocal thermal transfer.
    /// Account for actual endpoint elastic storage while the native return map
    /// uses frozen current-temperature stiffness. Energetic hardening stays fixed.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old: &ThermalCreepState<M::State>,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<ThermalCreepState<M::State>>, CreepError> {
        let old_thermal = self
            .potential
            .at_temperature(old.temperature_k)
            .map_err(thermal_error)?;
        let old_mechanical = self.mechanical_strain(old_strain, &old_thermal)?;
        let mechanical_strain = self.mechanical_strain(strain, &self.thermal)?;
        let native = self.mechanical.creep_increment(
            &old_mechanical,
            &old.mechanical,
            &mechanical_strain,
            tolerances,
        )?;
        let initial = self.thermal_state(old_strain, &old.mechanical, &old_thermal)?;
        let final_state = self.thermal_state(strain, &native.state, &self.thermal)?;
        let initial_stiffness = old_thermal.stiffness_mandel_pa();
        let final_stiffness = self.thermal.stiffness_mandel_pa();
        // Difference the stiffness first, not two large stored energies. This is
        // Psi(e_0,T_1)-Psi(e_0,T_0) with e_0 fixed, excluding free expansion.
        // Retain the native energetic-hardening contribution untouched.
        let coefficient_storage_change = 0.5
            * (0..6)
                .map(|i| {
                    initial.elastic_strain_mandel[i]
                        * (0..6)
                            .map(|j| {
                                (final_stiffness[i][j] - initial_stiffness[i][j])
                                    * initial.elastic_strain_mandel[j]
                            })
                            .sum::<f64>()
                })
                .sum::<f64>();
        let root_two = fs_math::det::sqrt(2.0);
        let stress_residual = final_state
            .stress_mandel_pa
            .iter()
            .zip(native.stress_pa)
            .enumerate()
            .map(|(i, (a, b))| (a - if i < 3 { b } else { root_two * b }).abs())
            .fold(0.0_f64, f64::max);
        if !stress_residual.is_finite() || stress_residual > tolerances.stress_pa {
            return Err(CreepError::StressResidual {
                residual_pa: stress_residual,
                tolerance_pa: tolerances.stress_pa,
            });
        }
        let selected = self
            .thermal
            .expansion
            .property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
            .ok_or(CreepError::Invalid("original expansion claim is missing"))?
            .answer()
            .receipt
            .selected;
        // Integrate THIS physical thermal interval through the original source
        // owner, avoiding subtraction of two large reference primitives.
        let delta_free = integrate_isotropic_thermal_expansion(
            self.potential.card(),
            &old_thermal.expansion.query(),
            &self.thermal.expansion.query(),
            MaterialPropertySelection::PinnedByProperty(vec![(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY.to_owned(),
                selected,
            )]),
        )
        .map_err(CreepError::Source)?
        .free_linear_strain();
        let mut energy = native.energy;
        energy.initial_stored_j_m3 -= coefficient_storage_change;
        energy.mechanical_work_j_m3 = contract(
            &native.stress_pa,
            &core::array::from_fn(|i| strain[i] - old_strain[i]),
        );
        energy.thermoelastic_storage_change_j_m3 = self.thermal.temperature_k()
            * final_state.entropy_j_k
            - old.temperature_k * initial.entropy_j_k;
        energy.reversible_heat_j_m3 = trace(&native.stress_pa) * delta_free
            - coefficient_storage_change
            - energy.thermoelastic_storage_change_j_m3;
        energy.residual_j_m3 = energy.mechanical_work_j_m3
            - (energy.final_stored_j_m3 - energy.initial_stored_j_m3)
            - energy.thermoelastic_storage_change_j_m3
            - energy.caloric_heat_j_m3()
            - energy.numerical_dissipation_j_m3;
        if [
            energy.mechanical_work_j_m3,
            energy.initial_stored_j_m3,
            energy.final_stored_j_m3,
            energy.thermoelastic_storage_change_j_m3,
            energy.reversible_heat_j_m3,
            energy.caloric_heat_j_m3(),
            energy.residual_j_m3,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || energy.residual_j_m3.abs() > tolerances.energy_j_m3
        {
            return Err(CreepError::EnergyResidual {
                residual_j_m3: energy.residual_j_m3,
                tolerance_j_m3: tolerances.energy_j_m3,
            });
        }
        Ok(CreepIncrement {
            strain: *strain,
            stress_pa: native.stress_pa,
            tangent_pa: native.tangent_pa,
            state: ThermalCreepState {
                mechanical: native.state,
                temperature_k: self.thermal.temperature_k(),
            },
            energy,
        })
    }
}

impl<M: CreepLaw> SmallStrainLaw for ThermalCreep<M>
where
    M::State: CreepHistoryStrain,
{
    type State = ThermalCreepState<M::State>;
    fn initial_state(&self) -> Self::State {
        ThermalCreepState {
            mechanical: self.mechanical.initial_state(),
            temperature_k: self.thermal.temperature_k(),
        }
    }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.trial_strain(strain, state)
            .map(|e| self.mechanical.stress(&e, &state.mechanical))
            .unwrap_or([f64::NAN; 6])
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.trial_strain(strain, state)
            .map(|e| self.mechanical.tangent(&e, &state.mechanical))
            .unwrap_or([[f64::NAN; 6]; 6])
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        let strain = self.trial_strain(strain, state).unwrap_or([f64::NAN; 6]);
        ThermalCreepState {
            mechanical: self.mechanical.update_state(&strain, &state.mechanical),
            temperature_k: if strain.iter().all(|x| x.is_finite()) {
                self.thermal.temperature_k()
            } else {
                f64::NAN
            },
        }
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.mechanical.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.mechanical.card()
    }
}
impl<M: CreepLaw> CreepLaw for ThermalCreep<M>
where
    M::State: CreepHistoryStrain,
{
    fn kinetics(&self) -> &NortonCreep {
        self.mechanical.kinetics()
    }
    fn validate_creep_history(
        &self,
        strain: &Voigt,
        state: &Self::State,
    ) -> Result<(), CreepError> {
        let thermal = self
            .potential
            .at_temperature(state.temperature_k)
            .map_err(thermal_error)?;
        let mechanical = self.mechanical_strain(strain, &thermal)?;
        self.mechanical
            .validate_creep_history(&mechanical, &state.mechanical)?;
        self.thermal_state(strain, &state.mechanical, &thermal)?;
        Ok(())
    }
    fn creep_increment(
        &self,
        old_strain: &Voigt,
        old: &Self::State,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<Self::State>, CreepError> {
        self.advance_with_energy(old_strain, old, strain, tolerances)
    }
}

fn thermal_error(error: ThermoelasticError) -> CreepError {
    match error {
        ThermoelasticError::Source(e) => CreepError::Source(e),
        ThermoelasticError::Material(e) => CreepError::Elastic(e),
        ThermoelasticError::Invalid(reason) => CreepError::Invalid(reason),
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CreepSourceError> {
    if cancelled() {
        Err(CreepSourceError::Cancelled)
    } else {
        Ok(())
    }
}
