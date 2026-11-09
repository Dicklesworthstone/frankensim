//! The original full elastic tensor as a native small-strain constitutive law.
//!
//! Tensor-Voigt strain/stress are related to Mandel coordinates by D and D^-1,
//! respectively, where D=diag(1,1,1,sqrt(2),sqrt(2),sqrt(2)). Therefore the
//! native tangent is D^-1 C_M D, not C_M: normal/shear couplings cannot be
//! copied unchanged. Work is epsilon:sigma = e_M dot s_M. See the notation
//! derivation in <https://csmbrannon.net/tag/mandel-notation/>.

use super::{ElasticAcousticState, invalid};
use crate::tensor::{from_mandel, mandel_stiffness_to_voigt, to_mandel};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
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
    fn admits_tensor_strain(&self, strain: &Voigt) -> bool {
        if let Some(limit) = self.frobenius_strain_limit() {
            // Use the source law's arithmetic on native tensor components.
            // An extra Mandel conversion can round an inclusive boundary
            // outward, so do not re-admit that converted strain afterward.
            crate::elastic::admits_elastic_strain(strain, limit)
        } else {
            self.admits_mandel_strain(&to_mandel(strain))
        }
    }

    // Linear action without strain admission: wave traces use unit amplitudes
    // and admit their physical superposition separately. The original K/G
    // modes also serve bulk mechanics; interfaces must not reconstruct them.
    pub(super) fn apply_stiffness_mandel(&self, strain: &Voigt) -> Voigt {
        if let Some((bulk, shear)) = self.isotropic_moduli_pa {
            let lambda = self.stiffness_mandel_pa[0][1];
            if lambda.abs() <= 0.5 * shear {
                // Near lambda=0 an uncoupled small component must not be
                // buried in subtraction of the spherical/deviatoric parts.
                // This preserves exact zero normal traction at grazing
                // incidence when the original tensor has lambda=0.
                let (trace, _) = crate::elastic::isotropic_strain_modes(strain);
                let coupling = lambda * trace;
                return std::array::from_fn(|i| {
                    (2.0 * shear).mul_add(strain[i], if i < 3 { coupling } else { 0.0 })
                });
            }
            return crate::elastic::isotropic_stress_from_moduli(bulk, 2.0 * shear, strain);
        }
        std::array::from_fn(|i| {
            self.stiffness_mandel_pa[i]
                .iter()
                .zip(strain)
                .map(|(c, e)| c * e)
                .sum()
        })
    }

    /// Evaluate the full anisotropic Hooke law and its recoverable energy.
    /// The original strain domain applies, including a native law's Frobenius
    /// norm bound. Normal/shear coupling and orientation are retained.
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
        self.strain_response_admitted(strain_mandel)
    }

    fn strain_response_admitted(
        &self,
        strain_mandel: Voigt,
    ) -> Result<ElasticStrainResponse, MaterialError> {
        let stress_mandel_pa = self.apply_stiffness_mandel(&strain_mandel);
        let stored_energy_density_j_m3 = if let Some((bulk, shear)) = self.isotropic_moduli_pa {
            let (trace, deviator) = crate::elastic::isotropic_strain_modes(&strain_mandel);
            // Positive modal storage retains a weak isochoric/hydrostatic
            // response instead of cancelling large, opposite stress-work terms.
            0.5 * (bulk * trace) * trace + deviator.iter().map(|e| (shear * e) * e).sum::<f64>()
        } else {
            0.5 * strain_mandel
                .iter()
                .zip(stress_mandel_pa)
                .map(|(e, s)| e * s)
                .sum::<f64>()
        };
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
        if !self.admits_tensor_strain(strain) {
            return [f64::NAN; 6];
        }
        match self.strain_response_admitted(to_mandel(strain)) {
            Ok(response) => from_mandel(&response.stress_mandel_pa),
            Err(_) => [f64::NAN; 6],
        }
    }

    fn tangent(&self, strain: &Voigt, (): &Self::State) -> Tangent6 {
        if !self.admits_tensor_strain(strain) {
            return [[f64::NAN; 6]; 6];
        }
        mandel_stiffness_to_voigt(&self.stiffness_mandel_pa)
    }

    fn update_state(&self, _strain: &Voigt, (): &Self::State) {}

    fn admits_strain(&self, strain: &Voigt) -> bool {
        self.admits_tensor_strain(strain)
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: Some(true),
            tangent_symmetric: true, // symmetric under tensor double contraction
            failure_envelope: "original elastic strain domain; no yield or damage law",
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
            // Exact for a native law's Frobenius ball, otherwise an inscribed
            // ball. Checked kernels retain the original domain shape.
            ValidityDomain::unconstrained().with(
                "strain-magnitude",
                0.0,
                self.maximum_absolute_mandel_strain,
            ),
            vec!["no finite strain, material evolution, yield or measured wave authority".into()],
            0.01,
        )
    }
}
