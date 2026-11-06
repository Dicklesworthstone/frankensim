//! Source-bound finite-strain elasticity with an explicitly Cauchy tensile limit.
//!
//! The existing hyperelastic law still owns energy, Piola stress and tangent.
//! This binding only pairs it with the caller's isotropic tensile criterion at
//! the same complete, fixed external material state. It supplies no damage,
//! fracture evolution, engineering-to-true-stress conversion or new source data.

use super::{IsotropicTensileLimit, invalid};
use crate::MaterialError;
use crate::hyper::Hyperelastic;
use crate::hyper::cards::{HyperelasticMaterial, ResolvedHyperelastic};
use crate::state_point::MaterialStatePointError;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::MaterialCard;

/// An admitted hyperelastic law and tensile limit from one immutable card.
///
/// The law, selected strength source and their original receipts remain intact.
/// A crossed limit is an observation; it never changes the constitutive law.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedTensileHyperelastic {
    elastic: ResolvedHyperelastic,
    limit: IsotropicTensileLimit,
    material_card_identity: ContentHash,
    identity: ContentHash,
}

impl ResolvedTensileHyperelastic {
    /// Bind an original model to an explicitly compatible Cauchy-strength limit.
    ///
    /// Calling this constructor asserts that the selected strength is suitable
    /// for comparison with true (Cauchy) stress. Engineering tensile strength,
    /// measured using original area, is not converted or silently interpreted
    /// as a finite-strain failure law. The isotropic cutoff is an explicit
    /// caller choice, even when the elastic model has directional structure.
    ///
    /// Both owners must come from this exact card and the same complete query,
    /// including coordinate quantities. Changing external conditions requires
    /// admitting both owners again before constructing a new binding.
    ///
    /// # Errors
    /// Refuses a foreign material card, an elastic model not present in that
    /// card, or any difference in the retained physical query.
    pub fn from_cauchy_limit(
        card: &MaterialCard,
        elastic: ResolvedHyperelastic,
        limit: IsotropicTensileLimit,
    ) -> Result<Self, MaterialStatePointError> {
        let material_card_identity = card.content_hash();
        if limit.resolved().card_identity() != material_card_identity {
            return Err(invalid(
                "finite-strain tensile strength requires the original material card",
            ));
        }
        let model_identity = elastic.model_card().content_hash();
        let model_index = card
            .models()
            .iter()
            .position(|model| {
                model == elastic.model_card() && model.content_hash() == model_identity
            })
            .ok_or_else(|| {
                invalid("finite-strain tensile elasticity requires an exact model from the material card")
            })?;
        if elastic.context() != &limit.resolved().query() {
            return Err(invalid(
                "finite-strain elasticity and tensile strength require the same complete typed query",
            ));
        }

        // The full card identity and exact member index bind the model without
        // copying or re-encoding its parameters. The bundle binds the complete
        // query and all source receipts; the selected property is also explicit.
        let mut hasher =
            DomainHasher::new("org.frankensim.fs-material.tensile-hyperelastic-cauchy.v1");
        hasher.update(material_card_identity.as_bytes());
        hasher.update(&(model_index as u64).to_le_bytes());
        hasher.update(limit.resolved().identity().as_bytes());
        hasher.update(&(limit.limit_property().len() as u64).to_le_bytes());
        hasher.update(limit.limit_property().as_bytes());
        hasher.update(&limit.limit_pa().to_bits().to_le_bytes());
        Ok(Self {
            elastic,
            limit,
            material_card_identity,
            identity: hasher.finalize(),
        })
    }

    /// Original source-bound energy, Piola stress and tangent owner.
    #[must_use]
    pub const fn elastic(&self) -> &ResolvedHyperelastic {
        &self.elastic
    }

    /// Selected true-stress tensile criterion and original source receipts.
    #[must_use]
    pub const fn limit(&self) -> &IsotropicTensileLimit {
        &self.limit
    }

    /// Immutable material card containing both the model and strength source.
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.material_card_identity
    }

    /// Exact model/strength/condition binding under the Cauchy interpretation.
    /// This content identity is not a fracture or safety certificate.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

impl HyperelasticMaterial for ResolvedTensileHyperelastic {
    fn hyperelastic(&self) -> &Hyperelastic {
        HyperelasticMaterial::hyperelastic(&self.elastic)
    }

    fn energy_density(&self, deformation: &[f64; 9]) -> f64 {
        HyperelasticMaterial::energy_density(&self.elastic, deformation)
    }

    fn energy_density_increment(&self, before: &[f64; 9], after: &[f64; 9]) -> f64 {
        HyperelasticMaterial::energy_density_increment(&self.elastic, before, after)
    }

    fn first_piola(&self, deformation: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        HyperelasticMaterial::first_piola(&self.elastic, deformation)
    }

    fn consistent_tangent(&self, deformation: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        HyperelasticMaterial::consistent_tangent(&self.elastic, deformation)
    }
}
