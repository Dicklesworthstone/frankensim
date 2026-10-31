//! Original source coefficients for simultaneous permanent creep and yielding.
//!
//! Norton kinetics, elastic storage and energetic isotropic hardening are fixed
//! over the complete supported temperature span. The independently sourced yield
//! offset may vary with temperature: it changes dissipation, not stored energy.
//! All nine original claims are selected together and pinned for continuation.
//! A monotonic post-yield tangent is converted only by the original J2 owner.
//! Nonzero translated yield surfaces require a noncoaxial coupled return map;
//! this source family therefore does not accept the mixed-hardening basis.

use std::sync::Arc;

use super::{NortonJ2, NortonJ2State};
use crate::creep::source::{
    CreepSourceError, CreepSourceFamily, NortonCreepSource, NortonSourceQuery, norton_requirements,
};
use crate::creep::{CreepError, CreepIncrement, CreepLaw, CreepTolerances, NortonCreep};
use crate::plastic::source::{HardeningPropertyBasis, ResolvedJ2Plasticity, j2_requirements};
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarPropertyRequirement,
    YIELD_STRESS_PROPERTY, resolve_material_state_point,
};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::ModelCard;
use fs_matdb::{ClaimSelection, MaterialCard};

/// The actual kinetic, elastic, yield and hardening source roles. Missing
/// hardening data never supplies an implicit perfectly plastic material.
pub fn norton_j2_requirements(
    basis: HardeningPropertyBasis,
) -> Result<Vec<ScalarPropertyRequirement>, CreepSourceError> {
    if basis == HardeningPropertyBasis::MixedPlasticStrainModuli {
        return Err(CreepSourceError::Invalid(
            "combined radial creep/plasticity requires isotropic hardening",
        ));
    }
    let mut requirements = norton_requirements()?;
    for requirement in j2_requirements(basis)? {
        if !requirements
            .iter()
            .any(|entry| entry.name() == requirement.name())
        {
            requirements.push(requirement);
        }
    }
    Ok(requirements)
}

#[derive(Debug)]
struct Binding<'a> {
    card: &'a MaterialCard,
    creep: NortonCreepSource<'a>,
    reference: ResolvedMaterialStatePoint,
    basis: HardeningPropertyBasis,
    strain_limit: f64,
}

/// Immutable source preparation shared across cells, retries and continuation.
/// Cloning this family shares its admitted tables and original claim receipts.
#[derive(Clone, Debug)]
pub struct NortonJ2Source<'a> {
    binding: Arc<Binding<'a>>,
}

impl<'a> NortonJ2Source<'a> {
    /// Admit a complete original temperature path. Every hardening knot is
    /// checked, so equal endpoints cannot hide changing recoverable storage.
    /// Yield may be piecewise linear, but must remain finite and positive.
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
            &norton_j2_requirements(basis)?,
            query.selection,
        )?;
        let selection = MaterialPropertySelection::PinnedByProperty(
            norton_requirements()?
                .iter()
                .map(|requirement| {
                    let property = reference
                        .property(requirement.name())
                        .expect("combined source contains every Norton requirement");
                    (
                        requirement.name().into(),
                        property.answer().receipt.selected,
                    )
                })
                .collect(),
        );
        let creep = NortonCreepSource::new(
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
        for name in [YIELD_STRESS_PROPERTY, basis.property_name()] {
            poll(&mut cancelled)?;
            let property = reference
                .property(name)
                .expect("combined source contains the selected plastic role");
            if !property.has_only_fixed_external_conditions(Some(query.temperature_axis), false) {
                return Err(CreepSourceError::Invalid(
                    "plastic coefficients cannot freeze evolving stress, strain, rate or frequency",
                ));
            }
            let selected = property.answer().receipt.selected;
            let claim = card
                .claims()
                .claim(selected)
                .ok_or(CreepSourceError::Invalid(
                    "selected plastic coefficient is absent",
                ))?;
            let line = card
                .claims()
                .query_line_typed(
                    &claim.key,
                    query.temperature_axis,
                    query.lower,
                    query.upper,
                    ClaimSelection::Pinned(selected),
                )
                .map_err(|error| CreepSourceError::Support {
                    property: name.into(),
                    reason: error.to_string(),
                })?;
            for (index, (_, answer)) in line.knots.iter().enumerate() {
                if index % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                let value = answer.evidence.value.value;
                if name == YIELD_STRESS_PROPERTY {
                    if !value.is_finite() || value <= 0.0 {
                        return Err(CreepSourceError::Invalid(
                            "yield stress must stay finite and positive throughout the source span",
                        ));
                    }
                } else if value != property.value_si() {
                    return Err(CreepSourceError::Invalid(
                        "energetic hardening must be constant over the complete thermal span",
                    ));
                }
            }
        }
        let plastic =
            ResolvedJ2Plasticity::try_from_material(reference.clone(), basis, query.strain_limit)?;
        NortonJ2::new(
            creep.prepared_law(creep.parameters().temperature_range_k[0], 1.0)?,
            plastic.law().clone(),
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            binding: Arc::new(Binding {
                card,
                creep,
                reference,
                basis,
                strain_limit: query.strain_limit,
            }),
        })
    }

    /// All original claims at the lower admitted temperature.
    pub fn reference(&self) -> &ResolvedMaterialStatePoint {
        &self.binding.reference
    }

    pub fn hardening_basis(&self) -> HardeningPropertyBasis {
        self.binding.basis
    }

    /// The original Norton source used by the shared caloric admission.
    pub fn creep_source(&self) -> &NortonCreepSource<'a> {
        &self.binding.creep
    }

    /// Resolve one real duration and temperature, retaining all selected claims.
    /// This does not evolve a specimen, material history or thermal clock.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<ResolvedNortonJ2, CreepSourceError> {
        let native = self.binding.creep.prepared_law(temperature_k, duration_s)?;
        let point = self
            .reference()
            .query()
            .with(self.binding.creep.temperature_axis(), temperature_k)
            .map_err(|error| CreepSourceError::Support {
                property: self.binding.creep.temperature_axis().into(),
                reason: error.to_string(),
            })?;
        let material = self.reference().at_query_point(self.binding.card, &point)?;
        let plastic = ResolvedJ2Plasticity::try_from_material(
            material.clone(),
            self.binding.basis,
            self.binding.strain_limit,
        )?;
        let law = NortonJ2::new(native, plastic.law().clone())?;
        Ok(ResolvedNortonJ2 {
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

/// One immutable native combined law and its complete original source bundle.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedNortonJ2 {
    material: Arc<ResolvedMaterialStatePoint>,
    law: NortonJ2,
    basis: HardeningPropertyBasis,
}

impl ResolvedNortonJ2 {
    pub const fn law(&self) -> &NortonJ2 {
        &self.law
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    pub const fn hardening_basis(&self) -> HardeningPropertyBasis {
        self.basis
    }
}

impl SmallStrainLaw for ResolvedNortonJ2 {
    type State = NortonJ2State;
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

impl CreepLaw for ResolvedNortonJ2 {
    fn kinetics(&self) -> &NortonCreep {
        self.law.creep()
    }
    fn validate_creep_history(
        &self,
        strain: &Voigt,
        state: &Self::State,
    ) -> Result<(), CreepError> {
        self.law.validate_history(strain, state)
    }
    fn creep_increment(
        &self,
        old_strain: &Voigt,
        old: &Self::State,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<Self::State>, CreepError> {
        self.law
            .advance_with_energy(old_strain, old, strain, tolerances)
    }
}

impl<'a> CreepSourceFamily<'a> for NortonJ2Source<'a> {
    type Material = ResolvedNortonJ2;
    fn creep_source(&self) -> &NortonCreepSource<'a> {
        NortonJ2Source::creep_source(self)
    }
    fn at_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<Self::Material, CreepSourceError> {
        NortonJ2Source::at_temperature(self, temperature_k, duration_s)
    }
}
