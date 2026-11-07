//! Source-defined transverse-shear Maxwell arms on a bonded circular section.
//!
//! Each layer contributes kappa*A_i*delta_G to the force/engineering-shear
//! relation. Original shear roles, correction, times and receipts are retained;
//! the receiving beam owns its spatial strain interpolation and memory states.

use super::{
    ConcentricShearSection, GeneralizedMaxwell, LayerRelaxationProperties, MaterialStatePointError,
    admit_bindings, invalid, resolve_layer_spectrum,
};

#[cfg(test)]
mod tests;

/// One original shear-modulus arm integrated over corrected layer area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransverseRelaxationBranch {
    pub layer_index: usize,
    /// kappa*A_i*delta_G [N], multiplying engineering shear.
    pub relaxing_rigidity_n: f64,
    pub relaxation_time_s: f64,
}

/// Fixed-state shear spectra on the original section, without a second history.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricPronyTransverseSection {
    shear: ConcentricShearSection,
    properties: Vec<LayerRelaxationProperties>,
    spectra: Vec<GeneralizedMaxwell>,
    branches: Vec<TransverseRelaxationBranch>,
    omega_band_rad_s: (f64, f64),
}

impl ConcentricShearSection {
    /// Resolve independent layer G(t) without changing equilibrium shear
    /// stiffness, mass, geometry or the source-defined shear correction.
    /// Explicit empty spectra declare an elastic channel; missing layer/source
    /// assignments refuse. The shared section resolver admits the exact G
    /// parents, positive times, nonnegative increments and frequency support.
    /// This does not infer shear relaxation from a Young-modulus spectrum.
    pub fn with_prony_transverse(
        self,
        properties: &[LayerRelaxationProperties],
    ) -> Result<ConcentricPronyTransverseSection, MaterialStatePointError> {
        let count = admit_bindings(self.elastic(), properties)?;
        let mut spectra = Vec::with_capacity(properties.len());
        let mut branches = Vec::with_capacity(count);
        let mut band = (0.0_f64, f64::MAX);
        let mut instantaneous = self.shear_stiffness_n();
        for (layer_index, binding) in properties.iter().enumerate() {
            let law = resolve_layer_spectrum(
                self.elastic(),
                Some(&self),
                layer_index,
                binding,
                &mut band,
            )?;
            let area = self.elastic().layer_integrals()[layer_index].area_m2();
            for &(delta, tau) in &law.terms {
                let rigidity = (delta * area) * self.shear_correction();
                instantaneous += rigidity;
                if !rigidity.is_finite()
                    || (delta > 0.0 && rigidity == 0.0)
                    || !instantaneous.is_finite()
                {
                    return Err(invalid("transverse relaxation rigidity is unrepresentable"));
                }
                branches.push(TransverseRelaxationBranch {
                    layer_index,
                    relaxing_rigidity_n: rigidity,
                    relaxation_time_s: tau,
                });
            }
            spectra.push(law);
        }
        if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0 {
            return Err(invalid(
                "section relaxation frequency bands do not intersect",
            ));
        }
        Ok(ConcentricPronyTransverseSection {
            shear: self,
            properties: properties.to_vec(),
            spectra,
            branches,
            omega_band_rad_s: band,
        })
    }
}

impl ConcentricPronyTransverseSection {
    /// Original equilibrium shear section and ordered material receipts.
    pub const fn shear(&self) -> &ConcentricShearSection {
        &self.shear
    }
    /// Original source roles in physical layer order.
    pub fn properties(&self) -> &[LayerRelaxationProperties] {
        &self.properties
    }
    /// Independent original shear spectra [Pa, s], one per layer.
    pub fn spectra(&self) -> &[GeneralizedMaxwell] {
        &self.spectra
    }
    /// Corrected-area arm rigidities [N], ordered by layer and declaration.
    pub fn branches(&self) -> &[TransverseRelaxationBranch] {
        &self.branches
    }
    /// Intersection of the original coefficient applicability intervals [rad/s].
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        self.omega_band_rad_s
    }
}
