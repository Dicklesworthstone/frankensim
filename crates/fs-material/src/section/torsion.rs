//! Source-bound torsional storage, rotary inertia and Kelvin-Voigt shear loss.
//!
//! Circular, perfectly bonded annuli share twist per length. The elastic owner
//! already supplies GJ; rotary inertia is sum(rho J) and viscous torque is
//! sum(eta_s J) times twist-rate per length, with J = 2I. Shear viscosity is an
//! EXPLICIT property, not inferred from axial/bending viscosity or a loss ratio.
//! These coefficients serve shafts, rods and strings without an instrument law.
//! Fixed material states and small-strain circular torsion only: no winding,
//! interlayer slip, warping of noncircular sections or rotational fluid loading.

use super::{ConcentricShearSection, MaterialStatePointError, positive, invalid};
use fs_qty::{Density, DynViscosity, QuantitySpec};
use super::relaxation::{LayerRelaxationProperties, LayerRelaxationBranch,
    SectionPronySpectrum, SectionRelaxationAxis};
use crate::visco::GeneralizedMaxwell;

/// Original source section with rotary mass and a separately selected shear law.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricTorsionalSection {
    shear: ConcentricShearSection,
    viscosity_property: String,
    rotary_inertia_kg_m: f64,
    viscous_stiffness_n_m2_s: f64,
    omega_band_rad_s: (f64, f64),
}

impl ConcentricShearSection {
    /// Bind torsional loss and inertia from the SAME ordered material bundles.
    ///
    /// Every layer must supply the named nonnegative shear viscosity [Pa s],
    /// with an explicit angular-frequency validity interval. Density, selected
    /// shear modulus and viscosity must be constants at the held typed
    /// environment; frequency samples cannot become a broadband law.
    /// Their applicability intervals are intersected, not extrapolated.
    /// Explicit zero viscosity is supported, missing loss data is not zero.
    /// The shear correction affects transverse shear only, never circular GJ.
    ///
    /// # Errors
    /// Missing/ambiguous/wrong-quantity coefficients, nonconstant laws, invalid
    /// or disjoint source bands, and unrepresentable section integrals refuse.
    pub fn with_torsional_viscosity(
        self, property_name: &str,
    ) -> Result<ConcentricTorsionalSection, MaterialStatePointError> {
        let mut inertia = 0.0;
        let mut viscous = 0.0;
        let mut band = (0.0_f64, f64::MAX);
        for (layer, geometry) in self.elastic().layers().iter().zip(self.elastic().layer_integrals()) {
            let rho = layer.material.property(crate::state_point::DENSITY_PROPERTY)
                .ok_or_else(|| invalid("torsion requires the original density"))?;
            let eta = layer.material.property(property_name)
                .ok_or_else(|| invalid("each torsional layer needs an explicit shear viscosity"))?;
            if rho.requirement().quantity() != QuantitySpec::dimensional(Density::DIMS)
                || eta.requirement().quantity() != QuantitySpec::dimensional(DynViscosity::DIMS)
                || !eta.value_si().is_finite() || eta.value_si() < 0.0 {
                return Err(invalid("torsion needs positive density and nonnegative dimension-only shear viscosity"));
            }
            let interval = eta.angular_frequency_band_rad_s()
                .map_err(|_| invalid("invalid shear-viscosity frequency semantics"))?
                .ok_or_else(|| invalid("torsional viscosity needs an explicit frequency validity band"))?;
            band = (band.0.max(interval.0), band.1.min(interval.1));
            for key in [crate::state_point::DENSITY_PROPERTY, self.shear_property(), property_name] {
                let property = layer.material.property(key)
                    .ok_or_else(|| invalid("torsion lost an original coefficient"))?;
                if !property.is_constant_at_fixed_environment() {
                    return Err(invalid("torsion requires fixed-state constant density, shear modulus and viscosity"));
                }
                if let Some((lo, hi)) = property.angular_frequency_band_rad_s()
                    .map_err(|_| invalid("invalid torsional coefficient frequency semantics"))? {
                    band = (band.0.max(lo), band.1.min(hi));
                }
            }
            let polar = positive(2.0 * geometry.second_moment_m4(), "polar area moment")?;
            let mass = positive(rho.value_si() * polar, "layer rotary inertia per length")?;
            inertia = positive(inertia + mass, "summed rotary inertia per length")?;
            let loss = eta.value_si() * polar;
            if !loss.is_finite() || (eta.value_si() > 0.0 && loss == 0.0) {
                return Err(invalid("torsional viscous stiffness is unrepresentable"));
            }
            viscous += loss;
            if !viscous.is_finite() { return Err(invalid("torsional loss sum overflow")); }
        }
        if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0 {
            return Err(invalid("torsional coefficient frequency bands do not intersect"));
        }
        Ok(ConcentricTorsionalSection { shear: self, viscosity_property: property_name.into(),
            rotary_inertia_kg_m: inertia, viscous_stiffness_n_m2_s: viscous, omega_band_rad_s: band })
    }
}

impl ConcentricTorsionalSection {
    /// Original shear, density, elastic section and every ordered source receipt.
    #[must_use]
    pub const fn shear(&self) -> &ConcentricShearSection { &self.shear }
    /// Explicit shear-viscosity selector, not a bending-loss alias.
    #[must_use]
    pub fn viscosity_property(&self) -> &str { &self.viscosity_property }
    /// Rotary inertia per unit physical length [kg m].
    #[must_use]
    pub const fn rotary_inertia_kg_m(&self) -> f64 { self.rotary_inertia_kg_m }
    /// Original circular Saint-Venant stiffness GJ [N m2].
    #[must_use]
    pub fn stiffness_n_m2(&self) -> f64 { self.shear.torsional_stiffness_n_m2() }
    /// Torque coefficient multiplying spatial twist-rate gradient [N m2 s].
    #[must_use]
    pub const fn viscous_stiffness_n_m2_s(&self) -> f64 { self.viscous_stiffness_n_m2_s }
    /// Intersection of the actual fixed-state material-law bands [rad/s].
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.omega_band_rad_s }
}

/// Frozen circular torsion with independent per-layer shear-relaxation histories.
/// The existing Kelvin dashpot remains a separately sourced parallel mechanism.
#[derive(Clone, Debug)]
pub struct ConcentricPronyTorsionalSection {
    torsional: ConcentricTorsionalSection,
    spectrum: SectionPronySpectrum,
}

impl ConcentricTorsionalSection {
    /// Bind each layer's shear spectrum using its actual polar area moment.
    ///
    /// The selected equilibrium modulus must match that layer's original G;
    /// a bending/Young modulus cannot silently replace it. A positive arm adds
    /// `delta_G * J` without changing equilibrium GJ or rotary inertia.
    /// Original time scales and source identities survive independently.
    /// Explicit empty spectra are elastic; omitted layers and source claims refuse.
    /// The existing shear viscosity supplies an independent parallel dashpot.
    ///
    /// This is fixed-state small-angle circular torsion for rods, shafts or
    /// strings. It does not supply axial creep, winding slip, noncircular
    /// warping, evolving temperature or a finite-strain constitutive model.
    ///
    /// # Errors
    /// The shared section-spectrum coefficient, frequency-support and arm-budget
    /// refusals apply. The original torsional material band must also intersect.
    pub fn with_prony_torsion(
        self,
        properties: &[LayerRelaxationProperties],
    ) -> Result<ConcentricPronyTorsionalSection, MaterialStatePointError> {
        let mut spectrum = SectionPronySpectrum::resolve(self.shear.elastic(),
            self.shear.shear_property(), self.stiffness_n_m2(),
            SectionRelaxationAxis::Torsion, properties)?;
        spectrum.omega_band_rad_s = (
            spectrum.omega_band_rad_s.0.max(self.omega_band_rad_s.0),
            spectrum.omega_band_rad_s.1.min(self.omega_band_rad_s.1),
        );
        if spectrum.omega_band_rad_s.0 > spectrum.omega_band_rad_s.1 {
            return Err(invalid("torsional relaxation and dashpot frequency bands do not intersect"));
        }
        Ok(ConcentricPronyTorsionalSection { torsional: self, spectrum })
    }
}

impl ConcentricPronyTorsionalSection {
    /// Original GJ, rotary inertia, independent Kelvin loss and all material receipts.
    #[must_use]
    pub const fn torsional(&self) -> &ConcentricTorsionalSection { &self.torsional }
    /// Source roles in physical core/sleeve order.
    #[must_use]
    pub fn properties(&self) -> &[LayerRelaxationProperties] { &self.spectrum.properties }
    /// Independent `delta_G * J` arms in layer/declaration order.
    #[must_use]
    pub fn branches(&self) -> &[LayerRelaxationBranch] { &self.spectrum.branches }
    /// Original shear-modulus spectra [Pa, s], one per layer.
    #[must_use]
    pub fn spectra(&self) -> &[GeneralizedMaxwell] { &self.spectrum.spectra }
    /// Intersection of all original elastic, viscous and relaxation supports.
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.spectrum.omega_band_rad_s }
}
