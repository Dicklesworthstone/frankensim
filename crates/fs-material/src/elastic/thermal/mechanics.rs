//! Native small-strain mechanics of the original source thermal potential.
//!
//! At fixed temperature, sigma=C(T)(epsilon-epsilon_free(T)). Preparation
//! retains a borrow of that complete source response and converts its tangent
//! once. Newton evaluations do not repeat thermal derivatives, source queries,
//! matrix inversions or allocations. The original `state` method still owns
//! free energy, entropy, internal energy and thermal derivatives.

use super::{ResolvedThermoelastic, ThermoelasticError};
use crate::tensor::{from_mandel, mandel_stiffness_to_voigt, to_mandel};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};

/// A fixed-temperature mechanical view of one original source response.
/// Both total and elastic Mandel strain retain the original component bounds.
/// Use per-element views for nonuniform temperature. A subsequent temperature
/// needs a newly resolved source response; this view never evolves heat or
/// accounts for thermal parameter work by itself.
#[derive(Clone, Debug)]
pub struct FrozenThermoelasticLaw<'a> {
    source: &'a ResolvedThermoelastic,
    stiffness_voigt_pa: Tangent6,
    free_strain_mandel: Voigt,
}

impl ResolvedThermoelastic {
    /// Prepare the same elastic potential for native plane, tetrahedral and
    /// hexahedral `SmallStrainLaw` consumers. In plane stress, the native local
    /// reduction solves all three transverse strains, including thermal shear.
    pub fn try_mechanical_law(&self) -> Result<FrozenThermoelasticLaw<'_>, ThermoelasticError> {
        let stiffness_voigt_pa = mandel_stiffness_to_voigt(&self.stiffness_mandel_pa());
        let free_strain_mandel = self.free_strain_mandel();
        if stiffness_voigt_pa
            .iter()
            .flatten()
            .chain(free_strain_mandel.iter())
            .any(|x| !x.is_finite())
        {
            return Err(ThermoelasticError::Invalid(
                "native thermal mechanical coefficients are unrepresentable",
            ));
        }
        Ok(FrozenThermoelasticLaw {
            source: self,
            stiffness_voigt_pa,
            free_strain_mandel,
        })
    }
}

impl FrozenThermoelasticLaw<'_> {
    /// Complete original response, including pinned elastic/expansion sources
    /// and its common thermodynamic `state(total_mandel_strain, volume)` owner.
    pub const fn source(&self) -> &ResolvedThermoelastic {
        self.source
    }

    fn elastic_strain(&self, total: &Voigt) -> Option<Voigt> {
        let total = to_mandel(total);
        let bound = self.source.maximum_absolute_mandel_strain;
        if total.iter().any(|e| !e.is_finite() || e.abs() > bound) {
            return None;
        }
        let elastic: Voigt =
            std::array::from_fn(|i| total[i] - self.free_strain_mandel[i]);
        if elastic.iter().any(|e| !e.is_finite() || e.abs() > bound) {
            return None;
        }
        Some(from_mandel(&elastic))
    }
}

impl SmallStrainLaw for FrozenThermoelasticLaw<'_> {
    type State = ();

    fn initial_state(&self) {}

    fn stress(&self, total: &Voigt, (): &Self::State) -> Voigt {
        let Some(elastic) = self.elastic_strain(total) else {
            return [f64::NAN; 6];
        };
        let stress: Voigt = std::array::from_fn(|i| {
            self.stiffness_voigt_pa[i]
                .iter()
                .zip(elastic)
                .map(|(c, e)| c * e)
                .sum()
        });
        let energy_density = 0.5 * crate::tensor::contract(&elastic, &stress);
        if stress.iter().any(|s| !s.is_finite())
            || !energy_density.is_finite()
            || energy_density < 0.0
        {
            [f64::NAN; 6]
        } else {
            stress
        }
    }

    fn tangent(&self, total: &Voigt, (): &Self::State) -> Tangent6 {
        if self.elastic_strain(total).is_some() {
            self.stiffness_voigt_pa
        } else {
            [[f64::NAN; 6]; 6]
        }
    }

    fn update_state(&self, _total: &Voigt, (): &Self::State) {}

    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: Some(true),
            tangent_symmetric: true,
            failure_envelope: "original total and elastic Mandel strain bounds; fixed temperature",
        }
    }

    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.frozen-source-thermoelastic",
            "0.1.0",
            Ambition::Solid,
            vec![
                "infinitesimal strain; original fixed-temperature thermal potential".into(),
                format!(
                    "elastic source {:?}; expansion source {:?}; temperature {} K",
                    self.source.elastic.identity(),
                    self.source.expansion.identity(),
                    self.source.temperature_k,
                ),
            ],
            ValidityDomain::unconstrained().with(
                "strain-magnitude",
                0.0,
                core::f64::consts::FRAC_1_SQRT_2 * self.source.maximum_absolute_mandel_strain,
            ),
            vec![
                "temperature evolution and reversible thermal work belong to the original thermal owner".into(),
            ],
            0.01,
        )
    }
}
