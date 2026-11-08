//! Source-bound PCM-composite conduction on the shared implicit heat owner.
//!
//! Calorics use fixed constituent MASS fractions; conductivity uses their actual
//! VOLUME fractions. Parallel/series layers are explicitly selected directional
//! homogenizations, not a universal foam or isotropic-composite mixing rule.
//! Original constituent density and conductivity changes enter the tangent even
//! on an isothermal latent plateau. Storage remains the original composite H.
//!
//! Geometry and reference mass belong to the transport owner. This does not
//! model expansion work, interfacial thermal resistance, melt convection,
//! percolation, or nonequilibrium constituent temperatures.
//!
//! Composite morphology and parallel limit: Barako et al., J. Appl. Phys. 124, 145103 (2018),
//! doi:10.1063/1.5031914, section II C.

use super::phase::{PhaseConductivityError, PreparedPhaseConductivity};
use crate::phase::composite::{PcmComposite, PcmCompositeResponse};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// Direction of heat flow relative to explicitly admitted constituent layers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PcmConductionMorphology {
    /// Layers parallel to heat flow: k = (1-phi) k_matrix + phi k_PCM.
    Parallel,
    /// Layers normal to heat flow: 1/k = (1-phi)/k_matrix + phi/k_PCM.
    Series,
}

/// One allocation-free evaluation retaining the original constituent state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PcmConductivitySample {
    pub caloric: PcmCompositeResponse,
    pub conductivity_w_m_k: f64,
    pub conductivity_slope_kg_w_m_k_j: f64,
}

/// Borrowed original calorics and prepared material conductivity sources.
/// No replacement enthalpy chart, source resampling or extra latent account.
#[derive(Clone, Debug)]
pub struct PreparedPcmConductivity<'a, 'm> {
    caloric: &'a PcmComposite<'m>,
    matrix: &'a PreparedPhaseConductivity<'m>,
    pcm: &'a PreparedPhaseConductivity<'m>,
    morphology: PcmConductionMorphology,
    identity: ContentHash,
}

impl<'a, 'm> PreparedPcmConductivity<'a, 'm> {
    /// Bind the exact original constituents over the full composite domain.
    ///
    /// Source laws must cover each constituent's complete used enthalpy range.
    /// Their own source/phase/pressure admission remains intact. Compatibility
    /// of the two materials' pressure and work conditions follows the original
    /// composite contract; composition need not be the same for both materials.
    /// The morphology is an explicit numerical model, not measured-source data.
    pub fn new(
        caloric: &'a PcmComposite<'m>,
        matrix: &'a PreparedPhaseConductivity<'m>,
        pcm: &'a PreparedPhaseConductivity<'m>,
        morphology: PcmConductionMorphology,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PhaseConductivityError> {
        checkpoint()?;
        if !std::ptr::eq(matrix.caloric_curve(), caloric.matrix_curve())
            || !std::ptr::eq(pcm.caloric_curve(), caloric.pcm_curve())
        {
            return Err(PhaseConductivityError::Invalid(
                "composite conductivity must retain both exact original caloric objects",
            ));
        }
        let knots = caloric.knots();
        for knot in [knots[0], knots[knots.len() - 1]] {
            checkpoint()?;
            matrix.state_at_specific_enthalpy(knot.matrix_specific_enthalpy_j_kg)?;
            pcm.state_at_specific_enthalpy(knot.pcm_specific_enthalpy_j_kg)?;
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.pcm-conductivity.v1");
        hash.update(caloric.identity().as_bytes());
        hash.update(matrix.identity().as_bytes());
        hash.update(pcm.identity().as_bytes());
        hash.update(&[match morphology {
            PcmConductionMorphology::Parallel => 0,
            PcmConductionMorphology::Series => 1,
        }]);
        let prepared = Self {
            caloric,
            matrix,
            pcm,
            morphology,
            identity: hash.finalize(),
        };
        for h in caloric.enthalpy_range_j_kg() {
            prepared.state_at_specific_enthalpy(h)?;
        }
        checkpoint()?;
        Ok(prepared)
    }

    pub const fn caloric_law(&self) -> &'a PcmComposite<'m> {
        self.caloric
    }
    pub const fn constituents(&self) -> [&'a PreparedPhaseConductivity<'m>; 2] {
        [self.matrix, self.pcm]
    }
    pub const fn morphology(&self) -> PcmConductionMorphology {
        self.morphology
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Analytic derivative includes conductivity and original density feedback.
    /// k may change while T is fixed during melting. Source knots use their
    /// original one-sided derivatives; no finite-difference query is performed.
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<PcmConductivitySample, PhaseConductivityError> {
        let caloric = self.caloric.response_at_specific_enthalpy(h)?;
        let matrix_h = caloric.state.matrix.specific_enthalpy_j_kg;
        let matrix = if h == self.caloric.enthalpy_range_j_kg()[1]
            && caloric.matrix_enthalpy_derivative != 0.0
        {
            // Retain both the admitted left conductivity source and its exact
            // Cp(T) conversion to dk/dh, including an interior source kink.
            self.matrix.state_at_specific_enthalpy_from_left(matrix_h)?
        } else {
            self.matrix.state_at_specific_enthalpy(matrix_h)?
        };
        let pcm = self
            .pcm
            .state_at_specific_enthalpy(caloric.state.pcm.specific_enthalpy_j_kg)?;
        let phi = caloric.state.pcm_volume_fraction;
        let dphi = caloric.pcm_volume_fraction_slope_kg_j;
        if phi <= 0.0 || phi >= 1.0 {
            return Err(PhaseConductivityError::Unrepresentable(
                "both original constituent volume fractions must remain representable",
            ));
        }
        let km = matrix.conductivity_w_m_k();
        let kp = pcm.conductivity_w_m_k();
        let dkm = matrix.conductivity_enthalpy_derivative_kg_w_m_k_j()
            * caloric.matrix_enthalpy_derivative;
        let dkp =
            pcm.conductivity_enthalpy_derivative_kg_w_m_k_j() * caloric.pcm_enthalpy_derivative;
        let (k, dk) = match self.morphology {
            PcmConductionMorphology::Parallel => (
                (1.0 - phi) * km + phi * kp,
                (1.0 - phi) * dkm + phi * dkp + (kp - km) * dphi,
            ),
            PcmConductionMorphology::Series if km == 0.0 && kp == 0.0 => (0.0, 0.0),
            PcmConductionMorphology::Series if km == 0.0 => (0.0, dkm / (1.0 - phi)),
            PcmConductionMorphology::Series if kp == 0.0 => (0.0, dkp / phi),
            PcmConductionMorphology::Series => {
                // Scale by the smaller conductivity, avoiding reciprocal
                // overflow for insulating but positive original constituents.
                let scale = km.min(kp);
                let k = scale / ((1.0 - phi) * (scale / km) + phi * (scale / kp));
                let (rm, rp) = (k / km, k / kp);
                let dk = (1.0 - phi) * rm * rm * dkm + phi * rp * rp * dkp + k * (rm - rp) * dphi;
                (k, dk)
            }
        };
        if ![dkm, dkp, k, dk].into_iter().all(f64::is_finite) || k < 0.0 || (k == 0.0 && dk != 0.0)
        {
            return Err(PhaseConductivityError::Unrepresentable(
                "composite conductivity or its original-state tangent is unrepresentable",
            ));
        }
        Ok(PcmConductivitySample {
            caloric,
            conductivity_w_m_k: k,
            conductivity_slope_kg_w_m_k_j: dk,
        })
    }
}

impl ConductivityLaw for PreparedPcmConductivity<'_, '_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.state_at_specific_enthalpy(h)
            .map(|s| s.conductivity_w_m_k)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}
impl ImplicitConductionLaw for PreparedPcmConductivity<'_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.caloric.enthalpy_range_j_kg()
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let s = self
            .state_at_specific_enthalpy(h)
            .map_err(|_| DepositionError::MaterialDomain)?;
        Ok(ImplicitConductionState {
            temperature_k: s.caloric.state.temperature_k,
            temperature_slope_kg_k_j: s.caloric.temperature_slope_kg_k_j,
            conductivity_w_m_k: s.conductivity_w_m_k,
            conductivity_slope_kg_w_m_k_j: s.conductivity_slope_kg_w_m_k_j,
        })
    }
    fn time_observables(&self, h: f64) -> Result<ImplicitConductionObservables, DepositionError> {
        let state = self
            .caloric
            .values_at_specific_enthalpy(h)
            .map_err(|_| DepositionError::MaterialDomain)?;
        Ok(ImplicitConductionObservables {
            temperature_k: state.temperature_k,
            phase_mass_fraction: Some(state.liquid_mass_fraction),
        })
    }
}
