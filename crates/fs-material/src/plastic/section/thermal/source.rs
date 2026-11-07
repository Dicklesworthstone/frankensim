//! Shared thermal-section inputs for native plasticity and permanent creep.
//!
//! This boundary carries original histories and native work accounts without
//! projecting creep into a plastic history. It prepares a law once per layer
//! and physical interval; the existing section and shell owners retain their
//! transverse condensation, quadrature, equilibrium and thermal integrators.

use crate::creep::source::{CreepSourceError, CreepSourceFamily};
use crate::creep::thermal::{FeedbackError as CreepFeedbackError, admit_zero_expansion_caloric};
use crate::creep::{CreepLaw, CreepTolerances};
use crate::phase::{EquilibriumPhaseState, ResolvedHeatCapacityCurve};
use crate::plastic::energy::EnergyTolerances;
use crate::plastic::hill::thermal::ThermalHillError;
use crate::plastic::thermal::family::{ThermoplasticFamilyError, ThermoplasticSourceFamily};
use crate::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use crate::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use crate::{SmallStrainLaw, Voigt};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};

/// One endpoint law at a calorically resolved temperature and physical duration.
#[derive(Clone, Debug)]
pub struct ThermalSectionStep<L> {
    pub temperature_k: f64,
    pub law: L,
}

/// Native history and per-reference-volume work, with no assumed J2 projection.
#[derive(Clone, Debug)]
pub struct ThermalSectionIncrement<S> {
    pub state: S,
    pub stress_pa: Voigt,
    /// Complete native physical loss; includes `creep_heat_j_m3` exactly once.
    pub physical_heat_j_m3: f64,
    pub creep_heat_j_m3: f64,
    pub reversible_heat_j_m3: f64,
    pub stored_change_j_m3: f64,
    pub coefficient_storage_change_j_m3: f64,
    pub thermoelastic_storage_change_j_m3: f64,
    pub numerical_loss_j_m3: f64,
    pub mechanical_work_j_m3: f64,
}

/// Preserve the native source/constitutive refusal at the shared section boundary.
#[derive(Debug)]
pub enum ThermalSectionSourceError {
    Plastic(ThermoplasticFamilyError),
    Hill(ThermalHillError),
    Creep(CreepFeedbackError),
    Transport(ThermalMaterialError),
    Invalid(&'static str),
}
impl core::fmt::Display for ThermalSectionSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Plastic(e) => e.fmt(f),
            Self::Hill(e) => e.fmt(f),
            Self::Creep(e) => e.fmt(f),
            Self::Transport(e) => e.fmt(f),
            Self::Invalid(what) => write!(f, "thermal section source: {what}"),
        }
    }
}
impl std::error::Error for ThermalSectionSourceError {}
impl From<ThermoplasticFamilyError> for ThermalSectionSourceError {
    fn from(error: ThermoplasticFamilyError) -> Self {
        Self::Plastic(error)
    }
}
impl From<ThermalHillError> for ThermalSectionSourceError {
    fn from(error: ThermalHillError) -> Self {
        Self::Hill(error)
    }
}
impl From<CreepFeedbackError> for ThermalSectionSourceError {
    fn from(error: CreepFeedbackError) -> Self {
        Self::Creep(error)
    }
}
impl From<CreepSourceError> for ThermalSectionSourceError {
    fn from(error: CreepSourceError) -> Self {
        Self::Creep(CreepFeedbackError::Source(error))
    }
}

/// Exact native source family consumed by layered thermal sections and shells.
/// Old caloric state/history and duration are frozen for every candidate.
pub trait ThermalSectionSourceFamily<'a>: Clone + core::fmt::Debug {
    type Material: SmallStrainLaw + Clone + core::fmt::Debug;

    fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve;
    fn reference_density_kg_m3(&self) -> f64;
    /// Borrow the original point no longer than either the family or its sources.
    fn thermal_material_point<'s>(&'s self) -> ThermalMaterialPoint<'s>
    where
        'a: 's;
    fn admit_section_enthalpy(
        &self,
        h: f64,
    ) -> Result<EquilibriumPhaseState, ThermalSectionSourceError>;
    fn resolve_section_step(
        &self,
        old_h: f64,
        new_h: f64,
        duration_s: f64,
    ) -> Result<ThermalSectionStep<Self::Material>, ThermalSectionSourceError>;
    /// Inspect the original accepted history at `old_h`, then audit exactly
    /// one native increment. The returned state retains every native memory.
    #[allow(clippy::too_many_arguments)]
    fn section_increment(
        &self,
        old_h: f64,
        new_h: f64,
        old_strain: &Voigt,
        old: &<Self::Material as SmallStrainLaw>::State,
        new_strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: f64,
    ) -> Result<
        ThermalSectionIncrement<<Self::Material as SmallStrainLaw>::State>,
        ThermalSectionSourceError,
    >;
    fn validate_section_history(
        &self,
        h: f64,
        strain: &Voigt,
        state: &<Self::Material as SmallStrainLaw>::State,
        stress_tolerance_pa: f64,
    ) -> Result<(), ThermalSectionSourceError>;

    /// Pin transport to the same original card/context and actual caloric state.
    /// Frequency, unowned strain and other evolving-coordinate fits refuse.
    fn pin_thermal_transport(
        &self,
        h: f64,
        key: &PropertyKey,
        selection: MaterialPropertySelection,
    ) -> Result<MaterialPropertySelection, ThermalSectionSourceError> {
        let state = self.admit_section_enthalpy(h)?;
        let thermal = self.thermal_material_point();
        let point = thermal
            .at(state)
            .map_err(ThermalSectionSourceError::Transport)?;
        let requirement =
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative)
                .map_err(|e| {
                    ThermalSectionSourceError::Transport(ThermalMaterialError::Material(e))
                })?;
        let material =
            resolve_material_state_point(thermal.material, &point, &[requirement], selection)
                .map_err(|e| {
                    ThermalSectionSourceError::Transport(ThermalMaterialError::Material(e))
                })?;
        let property = &material.properties()[0];
        thermal
            .admit_static(property)
            .map_err(ThermalSectionSourceError::Transport)?;
        if !property.has_only_fixed_external_conditions(Some(thermal.temperature_axis), false)
            || property
                .angular_frequency_band_rad_s()
                .map_err(|e| {
                    ThermalSectionSourceError::Transport(ThermalMaterialError::Material(e))
                })?
                .is_some()
        {
            return Err(ThermalSectionSourceError::Invalid(
                "instantaneous transport requires only owned temperature and fixed external conditions, without a frequency band",
            ));
        }
        Ok(MaterialPropertySelection::PinnedByProperty(vec![(
            key.name().to_owned(),
            property.answer().receipt.selected,
        )]))
    }
}

// Existing native J2, Voce and Chaboche source families preserve their exact
// source resolution and energy arithmetic. No second plastic return is added.
impl<'a, S: ThermoplasticSourceFamily<'a>> ThermalSectionSourceFamily<'a> for S {
    type Material = S::Material;
    fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve {
        self.thermal_source().caloric_source()
    }
    fn reference_density_kg_m3(&self) -> f64 {
        self.thermal_source().density_kg_m3()
    }
    fn thermal_material_point<'s>(&'s self) -> ThermalMaterialPoint<'s>
    where
        'a: 's,
    {
        self.thermal_source().thermal_material_point()
    }
    fn admit_section_enthalpy(
        &self,
        h: f64,
    ) -> Result<EquilibriumPhaseState, ThermalSectionSourceError> {
        self.thermal_source()
            .resolve_enthalpy(h)
            .map(|(state, _)| state)
            .map_err(|e| ThermalSectionSourceError::Plastic(e.into()))
    }
    fn resolve_section_step(
        &self,
        old_h: f64,
        new_h: f64,
        duration_s: f64,
    ) -> Result<ThermalSectionStep<Self::Material>, ThermalSectionSourceError> {
        self.admit_section_enthalpy(old_h)?;
        let step = ThermoplasticSourceFamily::resolve_step(self, old_h, new_h, duration_s)?;
        Ok(ThermalSectionStep {
            temperature_k: step.temperature_k,
            law: step.law,
        })
    }
    fn section_increment(
        &self,
        old_h: f64,
        new_h: f64,
        old_strain: &Voigt,
        old: &<Self::Material as SmallStrainLaw>::State,
        new_strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: f64,
    ) -> Result<
        ThermalSectionIncrement<<Self::Material as SmallStrainLaw>::State>,
        ThermalSectionSourceError,
    > {
        self.validate_accepted(old_h, old_strain, old, tolerances.yield_pa)?;
        let increment = self.checked_increment(
            old_h, new_h, old_strain, old, new_strain, tolerances, duration_s,
        )?;
        let thermal = increment.thermal;
        Ok(ThermalSectionIncrement {
            state: increment.state,
            stress_pa: thermal.mechanical.stress_pa,
            physical_heat_j_m3: thermal.mechanical.energy.physical_dissipation_j_m3,
            creep_heat_j_m3: 0.0,
            reversible_heat_j_m3: thermal.reversible_heat_j_m3,
            stored_change_j_m3: thermal.stored_change_j_m3(),
            coefficient_storage_change_j_m3: thermal.coefficient_storage_change_j_m3,
            thermoelastic_storage_change_j_m3: thermal.thermoelastic_storage_change_j_m3,
            numerical_loss_j_m3: thermal.mechanical.energy.numerical_dissipation_j_m3,
            mechanical_work_j_m3: thermal.mechanical_work_j_m3,
        })
    }
    fn validate_section_history(
        &self,
        h: f64,
        strain: &Voigt,
        state: &<Self::Material as SmallStrainLaw>::State,
        stress_tolerance_pa: f64,
    ) -> Result<(), ThermalSectionSourceError> {
        self.validate_accepted(h, strain, state, stress_tolerance_pa)
            .map_err(Into::into)
    }
    fn pin_thermal_transport(
        &self,
        h: f64,
        key: &PropertyKey,
        selection: MaterialPropertySelection,
    ) -> Result<MaterialPropertySelection, ThermalSectionSourceError> {
        self.thermal_source()
            .pin_thermal_transport(h, key, selection)
            .map_err(|e| ThermalSectionSourceError::Plastic(e.into()))
    }
}

/// Original creep family bound to one fully solid caloric source and reference
/// density. Supports standalone/joint, primary/secondary and explicitly sourced
/// expansion through the original `CreepSourceFamily`; no coefficient is inferred.
#[derive(Clone, Debug)]
pub struct CaloricCreepSource<'a, S: CreepSourceFamily<'a>> {
    card: &'a MaterialCard,
    source: S,
    caloric: &'a ResolvedHeatCapacityCurve,
    context: QueryPoint,
    density_kg_m3: f64,
}
impl<'a, S: CreepSourceFamily<'a>> CaloricCreepSource<'a, S> {
    pub fn new(
        card: &'a MaterialCard,
        source: S,
        caloric: &'a ResolvedHeatCapacityCurve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ThermalSectionSourceError> {
        let density =
            admit_zero_expansion_caloric(card, source.creep_source(), caloric, &mut cancelled)?;
        let mechanical_capacity = source.minimum_thermal_capacity_j_m3_k()?;
        for (index, pair) in caloric.curve().knots().windows(2).enumerate() {
            if cancelled() {
                return Err(CreepFeedbackError::Cancelled.into());
            }
            // The original caloric chart owns the Cp interpolation. Include
            // the expansion family's conservative mechanical capacity bound.
            let minimum_cp =
                if let Some([a, b]) = caloric.curve().sensible_heat_capacities_on_interval(index) {
                    a.min(b).next_down()
                } else {
                    let dh = (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
                        .next_down();
                    let dt = (pair[1].temperature_k - pair[0].temperature_k).next_up();
                    (dh / dt).next_down()
                };
            let caloric_capacity = (density * minimum_cp).next_down();
            let combined = (caloric_capacity + mechanical_capacity).next_down();
            if !mechanical_capacity.is_finite()
                || !caloric_capacity.is_finite()
                || caloric_capacity <= 0.0
                || !combined.is_finite()
                || combined <= 0.0
            {
                return Err(ThermalSectionSourceError::Invalid(
                    "positive combined caloric and mechanical capacity required over the original source span",
                ));
            }
        }
        let context = source.creep_source().reference().query();
        Ok(Self {
            card,
            source,
            caloric,
            context,
            density_kg_m3: density,
        })
    }
    pub fn creep_source_family(&self) -> &S {
        &self.source
    }
    fn phase(&self, h: f64) -> Result<EquilibriumPhaseState, ThermalSectionSourceError> {
        self.caloric
            .curve()
            .state_at_specific_enthalpy(h)
            .map_err(|e| CreepFeedbackError::Caloric(e).into())
    }
}

impl<'a, S: CreepSourceFamily<'a>> ThermalSectionSourceFamily<'a> for CaloricCreepSource<'a, S> {
    type Material = S::Material;
    fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve {
        self.caloric
    }
    fn reference_density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    fn thermal_material_point<'s>(&'s self) -> ThermalMaterialPoint<'s>
    where
        'a: 's,
    {
        ThermalMaterialPoint {
            material: self.card,
            context: &self.context,
            temperature_axis: self.source.creep_source().temperature_axis(),
            liquid_fraction_axis: None,
        }
    }
    fn admit_section_enthalpy(
        &self,
        h: f64,
    ) -> Result<EquilibriumPhaseState, ThermalSectionSourceError> {
        let state = self.phase(h)?;
        self.source.at_temperature(state.temperature_k(), 1.0)?;
        Ok(state)
    }
    fn resolve_section_step(
        &self,
        old_h: f64,
        new_h: f64,
        duration_s: f64,
    ) -> Result<ThermalSectionStep<Self::Material>, ThermalSectionSourceError> {
        self.admit_section_enthalpy(old_h)?;
        let state = self.phase(new_h)?;
        let law = self
            .source
            .at_temperature(state.temperature_k(), duration_s)?;
        Ok(ThermalSectionStep {
            temperature_k: state.temperature_k(),
            law,
        })
    }
    fn section_increment(
        &self,
        old_h: f64,
        new_h: f64,
        old_strain: &Voigt,
        old: &<Self::Material as SmallStrainLaw>::State,
        new_strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: f64,
    ) -> Result<
        ThermalSectionIncrement<<Self::Material as SmallStrainLaw>::State>,
        ThermalSectionSourceError,
    > {
        self.validate_section_history(old_h, old_strain, old, tolerances.yield_pa)?;
        let phase = self.phase(new_h)?;
        let law = self
            .source
            .at_temperature(phase.temperature_k(), duration_s)?;
        let value = law
            .creep_increment(
                old_strain,
                old,
                new_strain,
                CreepTolerances {
                    stress_pa: tolerances.yield_pa,
                    energy_j_m3: tolerances.energy_j_m3,
                },
            )
            .map_err(CreepFeedbackError::Mechanical)?;
        let energy = value.energy;
        Ok(ThermalSectionIncrement {
            state: value.state,
            stress_pa: value.stress_pa,
            physical_heat_j_m3: energy.physical_dissipation_j_m3,
            creep_heat_j_m3: energy.creep_dissipation_j_m3,
            reversible_heat_j_m3: energy.reversible_heat_j_m3,
            stored_change_j_m3: energy.final_stored_j_m3 - energy.initial_stored_j_m3,
            // Native Norton families admit fixed elastic/hardening storage.
            coefficient_storage_change_j_m3: 0.0,
            thermoelastic_storage_change_j_m3: energy.thermoelastic_storage_change_j_m3,
            numerical_loss_j_m3: energy.numerical_dissipation_j_m3,
            mechanical_work_j_m3: energy.mechanical_work_j_m3,
        })
    }
    fn validate_section_history(
        &self,
        h: f64,
        strain: &Voigt,
        state: &<Self::Material as SmallStrainLaw>::State,
        stress_tolerance_pa: f64,
    ) -> Result<(), ThermalSectionSourceError> {
        if !stress_tolerance_pa.is_finite() || stress_tolerance_pa <= 0.0 {
            return Err(ThermalSectionSourceError::Invalid(
                "positive finite stress allowance required",
            ));
        }
        let phase = self.phase(h)?;
        self.source
            .validate_accepted_temperature(phase.temperature_k(), strain, state)?;
        Ok(())
    }
}
