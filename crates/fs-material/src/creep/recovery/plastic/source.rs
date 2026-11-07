//! Source-pinned recovery, permanent creep and isotropic plasticity.
//!
//! Fifteen original claims bind the complete simultaneous law. The original
//! recovery source owns its thirteen creep/resistance coefficients; the shared
//! joint source admission owns temperature-dependent yield and fixed energetic
//! hardening. No source refresh changes specimen history or advances time.

use std::sync::Arc;

use super::{RecoveringNortonJ2, RecoveringNortonJ2State};
use crate::creep::plastic::NortonJ2;
use crate::creep::plastic::source::{admit_plastic_span, combined_requirements};
use crate::creep::recovery::CreepRecoveryParameters;
use crate::creep::recovery::source::{RecoveringNortonSource, recovering_norton_requirements};
use crate::creep::source::{
    CreepSourceError, CreepSourceFamily, NortonCreepSource, NortonSourceQuery,
};
use crate::creep::{CreepError, CreepIncrement, CreepLaw, CreepTolerances, NortonCreep};
use crate::plastic::source::{HardeningPropertyBasis, ResolvedJ2Plasticity};
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::ModelCard;
use fs_matdb::MaterialCard;

/// All primary flow, resistance recovery, elastic, yield and hardening claims.
/// Mixed kinematic hardening requires a noncoaxial return and is refused.
pub fn recovering_norton_j2_requirements(
    basis: HardeningPropertyBasis,
) -> Result<Vec<ScalarPropertyRequirement>, CreepSourceError> {
    combined_requirements(basis, recovering_norton_requirements()?)
}

#[derive(Debug)]
struct Binding<'a> {
    card: &'a MaterialCard,
    recovery: RecoveringNortonSource<'a>,
    reference: ResolvedMaterialStatePoint,
    basis: HardeningPropertyBasis,
    strain_limit: f64,
}

/// Immutable preparation, shared across cells, intervals and thermal trials.
#[derive(Clone, Debug)]
pub struct RecoveringNortonJ2Source<'a> {
    binding: Arc<Binding<'a>>,
}

impl<'a> RecoveringNortonJ2Source<'a> {
    /// Select all fifteen original claims together. Each supported thermal knot
    /// is checked by the existing recovery and energetic-plasticity owners.
    /// Yield may vary with temperature; reference kinetics, E, nu and H stay fixed.
    pub fn new(
        card: &'a MaterialCard,
        query: NortonSourceQuery<'_>,
        basis: HardeningPropertyBasis,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CreepSourceError> {
        poll(&mut cancelled)?;
        let reference = resolve_material_state_point(
            card,
            query.lower,
            &recovering_norton_j2_requirements(basis)?,
            query.selection,
        )?;
        let selection = MaterialPropertySelection::PinnedByProperty(
            recovering_norton_requirements()?
                .into_iter()
                .map(|requirement| {
                    let property = reference
                        .property(requirement.name())
                        .expect("joint source contains every recovery coefficient");
                    (
                        requirement.name().into(),
                        property.answer().receipt.selected,
                    )
                })
                .collect(),
        );
        let recovery = RecoveringNortonSource::new(
            card,
            NortonSourceQuery {
                lower: query.lower,
                upper: query.upper,
                temperature_axis: query.temperature_axis,
                selection,
                strain_limit: query.strain_limit,
            },
            &mut cancelled,
        )?;
        admit_plastic_span(
            card,
            &reference,
            basis,
            query.lower,
            query.upper,
            query.temperature_axis,
            &mut cancelled,
        )?;
        let plastic =
            ResolvedJ2Plasticity::try_from_material(reference.clone(), basis, query.strain_limit)?;
        let native = recovery.creep_source().prepared_law(
            recovery.creep_source().parameters().temperature_range_k[0],
            1.0,
        )?;
        RecoveringNortonJ2::new(
            NortonJ2::new(native, plastic.law().clone())?,
            *recovery.recovery_parameters(),
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            binding: Arc::new(Binding {
                card,
                recovery,
                reference,
                basis,
                strain_limit: query.strain_limit,
            }),
        })
    }

    /// Complete original material/claim bundle at the lower source temperature.
    pub fn reference(&self) -> &ResolvedMaterialStatePoint {
        &self.binding.reference
    }

    pub fn creep_source(&self) -> &NortonCreepSource<'a> {
        self.binding.recovery.creep_source()
    }

    pub fn recovery_source(&self) -> &RecoveringNortonSource<'a> {
        &self.binding.recovery
    }

    pub fn recovery_parameters(&self) -> &CreepRecoveryParameters {
        self.binding.recovery.recovery_parameters()
    }

    pub fn hardening_basis(&self) -> HardeningPropertyBasis {
        self.binding.basis
    }

    /// Resolve one physical interval and temperature. The complete pinned bundle
    /// refreshes once; subsequent native trials do not query material tables.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<ResolvedRecoveringNortonJ2, CreepSourceError> {
        let native = self
            .creep_source()
            .prepared_law(temperature_k, duration_s)?;
        let point = self
            .reference()
            .query()
            .with(self.creep_source().temperature_axis(), temperature_k)
            .map_err(|error| CreepSourceError::Support {
                property: self.creep_source().temperature_axis().into(),
                reason: error.to_string(),
            })?;
        let material = self.reference().at_query_point(self.binding.card, &point)?;
        let plastic = ResolvedJ2Plasticity::try_from_material(
            material.clone(),
            self.binding.basis,
            self.binding.strain_limit,
        )?;
        let law = RecoveringNortonJ2::new(
            NortonJ2::new(native, plastic.law().clone())?,
            *self.recovery_parameters(),
        )?;
        Ok(ResolvedRecoveringNortonJ2 {
            material: Arc::new(material),
            law,
            basis: self.binding.basis,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CreepSourceError> {
    if cancelled() {
        Err(CreepSourceError::Cancelled)
    } else {
        Ok(())
    }
}

/// One original simultaneous law, separate histories and complete source receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedRecoveringNortonJ2 {
    material: Arc<ResolvedMaterialStatePoint>,
    law: RecoveringNortonJ2,
    basis: HardeningPropertyBasis,
}

impl ResolvedRecoveringNortonJ2 {
    pub const fn law(&self) -> &RecoveringNortonJ2 {
        &self.law
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    pub const fn hardening_basis(&self) -> HardeningPropertyBasis {
        self.basis
    }
}

impl SmallStrainLaw for ResolvedRecoveringNortonJ2 {
    type State = RecoveringNortonJ2State;
    fn initial_state(&self) -> Self::State {
        self.law.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.law.stress(strain, state)
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.law.tangent(strain, state)
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.law.update_state(strain, state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.law.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.law.card()
    }
}

impl CreepLaw for ResolvedRecoveringNortonJ2 {
    fn kinetics(&self) -> &NortonCreep {
        self.law.kinetics()
    }
    fn validate_creep_history(
        &self,
        strain: &Voigt,
        state: &Self::State,
    ) -> Result<(), CreepError> {
        self.law.validate_creep_history(strain, state)
    }
    fn creep_increment(
        &self,
        old_strain: &Voigt,
        old: &Self::State,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<Self::State>, CreepError> {
        self.law
            .creep_increment(old_strain, old, strain, tolerances)
    }
}

impl<'a> CreepSourceFamily<'a> for RecoveringNortonJ2Source<'a> {
    type Material = ResolvedRecoveringNortonJ2;
    fn creep_source(&self) -> &NortonCreepSource<'a> {
        RecoveringNortonJ2Source::creep_source(self)
    }
    fn at_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<Self::Material, CreepSourceError> {
        RecoveringNortonJ2Source::at_temperature(self, temperature_k, duration_s)
    }
}
