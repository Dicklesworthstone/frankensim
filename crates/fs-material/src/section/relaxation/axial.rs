//! Source-bound axial Maxwell arms for perfectly bonded circular members.
//!
//! Common axial strain gives equilibrium rigidity `sum(E_eq A)` and each
//! original material arm contributes `delta_E A` [N]. Area, rather than the
//! bending or polar area moment, is the work-conjugate section integral.
//! Original layer identities, source times and coefficient validity survive.
//! The receiving structural owner supplies strain histories and their evolution;
//! this immutable section does not invent an initial material state.

use super::{
    ConcentricElasticSection, GeneralizedMaxwell, LayerRelaxationProperties,
    MaterialStatePointError, admit_bindings, invalid, resolve_layer_spectrum,
};

/// One original Young-modulus arm integrated over its actual annular area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxialRelaxationBranch {
    /// Index of the original spatial layer; source order is preserved.
    pub layer_index: usize,
    /// Excess axial rigidity `delta_E * A` [N], multiplying branch elastic strain.
    /// This is additional to the section's EQUILIBRIUM axial rigidity.
    pub axial_stiffness_n: f64,
    /// Original positive source time [s], without spectrum averaging.
    pub relaxation_time_s: f64,
}

/// Original equilibrium section and passive axial material spectra.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricPronyAxialSection {
    elastic: ConcentricElasticSection,
    properties: Vec<LayerRelaxationProperties>,
    spectra: Vec<GeneralizedMaxwell>,
    branches: Vec<AxialRelaxationBranch>,
    omega_band_rad_s: (f64, f64),
}

impl ConcentricElasticSection {
    /// Bind each layer's independent Young-modulus Maxwell spectrum to axial
    /// extension, retaining equilibrium EA, inertia, geometry and source receipts.
    ///
    /// Every layer supplies an explicit binding. Its equilibrium modulus must
    /// equal the original elastic Young modulus; branch coefficients are held
    /// at the same fixed environment. An explicitly empty spectrum is elastic,
    /// while missing data refuse. Source bands intersect through the existing
    /// shared section resolver, and at most 128 total source arms are admitted.
    /// A rod of reference length L uses each `axial_stiffness_n / L` as the
    /// spring stiffness conjugate to its branch elastic extension.
    ///
    /// # Errors
    /// Incomplete or invalid source assignments, incompatible quantities or
    /// equilibrium moduli, unsupported coefficient dependence, negative modulus
    /// increments, nonpositive times, disjoint bands, or unrepresentable rigidity.
    pub fn with_prony_axial(
        self,
        properties: &[LayerRelaxationProperties],
    ) -> Result<ConcentricPronyAxialSection, MaterialStatePointError> {
        let count = admit_bindings(&self, properties)?;
        let mut spectra = Vec::with_capacity(properties.len());
        let mut branches = Vec::with_capacity(count);
        let mut band = (0.0_f64, f64::MAX);
        let mut instantaneous = self.axial_stiffness_n();
        for (layer_index, binding) in properties.iter().enumerate() {
            let law = resolve_layer_spectrum(&self, None, layer_index, binding, &mut band)?;
            let area = self.layer_integrals()[layer_index].area_m2();
            for &(delta, tau) in &law.terms {
                let rigidity = delta * area;
                instantaneous += rigidity;
                if !rigidity.is_finite()
                    || (delta > 0.0 && rigidity == 0.0)
                    || !instantaneous.is_finite()
                {
                    return Err(invalid("axial relaxation rigidity is unrepresentable"));
                }
                branches.push(AxialRelaxationBranch {
                    layer_index,
                    axial_stiffness_n: rigidity,
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
        Ok(ConcentricPronyAxialSection {
            elastic: self,
            properties: properties.to_vec(),
            spectra,
            branches,
            omega_band_rad_s: band,
        })
    }
}

impl ConcentricPronyAxialSection {
    /// Original equilibrium elasticity, geometry, mass and every material receipt.
    pub const fn elastic(&self) -> &ConcentricElasticSection {
        &self.elastic
    }
    /// Original property selections in spatial layer order.
    pub fn properties(&self) -> &[LayerRelaxationProperties] {
        &self.properties
    }
    /// Original independent Young spectra [Pa, s], one per layer.
    pub fn spectra(&self) -> &[GeneralizedMaxwell] {
        &self.spectra
    }
    /// Area-integrated excess rigidities, ordered by layer and source declaration.
    pub fn branches(&self) -> &[AxialRelaxationBranch] {
        &self.branches
    }
    /// Intersection of original coefficient applicability intervals [rad/s].
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        self.omega_band_rad_s
    }
}
