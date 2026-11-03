//! Competing dissolved metal pools coupled to finite gold cyanidation.
//!
//! Each supplied metal family shares the original HCN/CN acid family and its
//! source reference. Every rate stage re-solves finite ligand conservation,
//! metal distributions and electroneutrality. Binding is reversible equilibrium;
//! it is not irreversible cyanide destruction. Metal totals remain fixed while
//! cyanide is transferred into newly formed Au(CN)2-.
//!
//! A Cu(I) family must be supplied as Cu(I). Neither the presence of Cu(II)
//! source data nor an ore copper assay infers reduction, mineral dissolution,
//! accessible copper, or a source-supported reaction rate.

use super::*;
use fs_thermochem::stoichiometry::thermodynamics::acid_base::complexation::{
    IdealLigandSystem, LigandSpeciationOptions, MetalLigandFamily,
};

/// Explicit, finite dissolved metals competing for the supplied cyanide.
#[derive(Clone, Debug)]
pub struct CompetingMetals {
    /// Source-backed association chains; identities and sources must agree with
    /// the original cyanidation chemistry. At least one family is required.
    pub families: Vec<MetalLigandFamily>,
    /// Every metal family exactly once, including explicitly zero totals.
    /// These are dissolved metal amounts, not whole-ore assays.
    pub totals: Vec<(String, Amount)>,
    /// Bounded simultaneous ligand/pH solve controls at every kinetic stage.
    pub numerics: LigandSpeciationOptions,
}

#[derive(Clone, Debug)]
pub(super) struct BoundCompetition {
    pub model: IdealLigandSystem,
    pub totals: Vec<(String, Amount)>,
    pub numerics: LigandSpeciationOptions,
}

impl GoldCyanidationBatch {
    /// Initialize a bath with finite competing dissolved metal pools.
    ///
    /// `charge.uncomplexed_cyanide` supplies the ligand total before competition;
    /// `state().uncomplexed_cyanide` reports its actual HCN+CN remainder after
    /// binding. Explicit metal counterions belong in `charge.spectators`.
    /// Initial speciation and every kinetic source coordinate include metal
    /// competition from the outset. The same coupling survives attachment to
    /// the finite oxygen-transfer adapter.
    pub fn new_with_competing_metals(
        model: GoldCyanidation,
        charge: CyanidationCharge,
        conditions: QueryPoint,
        metals: CompetingMetals,
    ) -> Result<Self> {
        require(
            !metals.families.is_empty(),
            "competing metals require at least one declared association family",
        )?;
        require(
            metals.families.iter().all(|family| {
                family.species().iter().all(|species| {
                    !model
                        .reaction()
                        .terms()
                        .iter()
                        .any(|term| &term.species == species)
                })
            }),
            "competing metal members cannot duplicate the original leach reactant or product stocks",
        )?;
        let chemistry = IdealLigandSystem::new(
            model.chemistry().clone(),
            model.cyanide_family().to_owned(),
            metals.families,
        )?;
        let speciation = metals.numerics.aqueous;
        Self::new_inner(
            model,
            charge,
            conditions,
            speciation,
            Some(BoundCompetition {
                model: chemistry,
                totals: metals.totals,
                numerics: metals.numerics,
            }),
        )
    }
}
