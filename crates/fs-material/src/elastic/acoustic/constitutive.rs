//! The original full elastic tensor as a native small-strain constitutive law.
//!
//! Tensor-Voigt strain/stress are related to Mandel coordinates by D and D^-1,
//! respectively, where D=diag(1,1,1,sqrt(2),sqrt(2),sqrt(2)). Therefore the
//! native tangent is D^-1 C_M D, not C_M: normal/shear couplings cannot be
//! copied unchanged. Work is epsilon:sigma = e_M dot s_M. See the notation
//! derivation in <https://csmbrannon.net/tag/mandel-notation/>.

use super::{ElasticAcousticState, invalid};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use crate::tensor::{from_mandel, mandel_stiffness_to_voigt, to_mandel};
use fs_blake3::ContentHash;
use fs_evidence::{Ambition, ModelCard, ValidityDomain};

/// Checked reversible response of the same tensor used by mechanics and waves.
#[derive(Clone, Debug, PartialEq)]
pub struct ElasticStrainResponse {
    /// Original frozen material identity, without refitting or source selection.
    pub material_state_identity: ContentHash,
    /// Total strain in `[xx, yy, zz, sqrt(2)xy, sqrt(2)yz, sqrt(2)zx]` order.
    pub strain_mandel: Voigt,
    /// Conjugate stress in the same Mandel order [Pa].
    pub stress_mandel_pa: Voigt,
    /// Recoverable energy per reference volume, `epsilon:sigma/2` [J/m^3].
    pub stored_energy_density_j_m3: f64,
}

impl ElasticAcousticState {
    fn admits_mandel_strain(&self, strain: &Voigt) -> bool {
        strain
            .iter()
            .all(|v| v.is_finite() && v.abs() <= self.maximum_absolute_mandel_strain)
    }

    /// Evaluate the full anisotropic Hooke law and its recoverable energy.
    /// The original componentwise Mandel strain ceiling applies; neither
    /// normal/shear coupling nor orientation is replaced by isotropic data.
    /// No allocation, factorization or source lookup occurs during evaluation.
    pub fn strain_response(
        &self,
        strain_mandel: Voigt,
    ) -> Result<ElasticStrainResponse, MaterialError> {
        if !self.admits_mandel_strain(&strain_mandel) {
            return Err(invalid(
                "elastic response exceeds the original Mandel strain domain",
            ));
        }
        let stress_mandel_pa: Voigt = std::array::from_fn(|i| {
            self.stiffness_mandel_pa[i]
                .iter()
                .zip(strain_mandel)
                .map(|(c, e)| c * e)
                .sum()
        });
        let stored_energy_density_j_m3 = 0.5
            * strain_mandel
                .iter()
                .zip(stress_mandel_pa)
                .map(|(e, s)| e * s)
                .sum::<f64>();
        if stress_mandel_pa.iter().any(|s| !s.is_finite())
            || !stored_energy_density_j_m3.is_finite()
            || stored_energy_density_j_m3 < 0.0
        {
            return Err(invalid(
                "elastic stress or stored energy is unrepresentable",
            ));
        }
        Ok(ElasticStrainResponse {
            material_state_identity: self.material_state_identity,
            strain_mandel,
            stress_mandel_pa,
            stored_energy_density_j_m3,
        })
    }
}

/// A frozen full tensor can directly drive the existing plane-stress/strain,
/// tetrahedral and hexahedral `SmallStrainLaw` consumers. Their constrained
/// reductions use this exact tangent and retain all transverse couplings.
///
/// The trait's infallible kernels return nonfinite values outside the original
/// source strain domain, so checked solvers refuse instead of extending it.
/// Use `strain_response` for a structured error. This state has no plasticity,
/// thermal eigenstrain, history, dissipation or finite-rotation constitutive law.
impl SmallStrainLaw for ElasticAcousticState {
    type State = ();

    fn initial_state(&self) {}

    fn stress(&self, strain: &Voigt, (): &Self::State) -> Voigt {
        match self.strain_response(to_mandel(strain)) {
            Ok(response) => from_mandel(&response.stress_mandel_pa),
            Err(_) => [f64::NAN; 6],
        }
    }

    fn tangent(&self, strain: &Voigt, (): &Self::State) -> Tangent6 {
        if !self.admits_mandel_strain(&to_mandel(strain)) {
            return [[f64::NAN; 6]; 6];
        }
        mandel_stiffness_to_voigt(&self.stiffness_mandel_pa)
    }

    fn update_state(&self, _strain: &Voigt, (): &Self::State) {}

    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: Some(true),
            tangent_symmetric: true, // symmetric under tensor double contraction
            failure_envelope: "original Mandel strain ceiling; no yield or damage law",
        }
    }

    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.full-tensor-elastic",
            "0.1.0",
            Ambition::Solid,
            vec![
                "frozen symmetric positive-definite elastic tensor; infinitesimal strain".into(),
                format!("original material state {:?}", self.material_state_identity),
            ],
            // The generic evidence helper uses max absolute tensor-Voigt
            // strain. This conservative scalar ceiling also bounds Mandel
            // shear; checked kernels apply the exact componentwise domain.
            ValidityDomain::unconstrained().with(
                "strain-magnitude",
                0.0,
                core::f64::consts::FRAC_1_SQRT_2 * self.maximum_absolute_mandel_strain,
            ),
            vec!["no finite strain, material evolution, yield or measured wave authority".into()],
            0.01,
        )
    }
}
