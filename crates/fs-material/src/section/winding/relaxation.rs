//! Linear viscoelastic correspondence on the ORIGINAL close-coiled spring.
//!
//! Its axial rigidity is a fixed positive geometry factor times wire shear
//! modulus. Each declared Maxwell increment therefore contributes that SAME
//! factor times Delta G, with the original relaxation time. This adds no
//! inter-turn friction, wrapping prestress, slip law or winding inertia.

use super::*;
use crate::section::relaxation::{LayerRelaxationProperties, MAX_SECTION_RELAXATION_BRANCHES};
use crate::visco::GeneralizedMaxwell;
use fs_qty::Time;

/// One original shear arm expressed in work-conjugate axial rigidity [N].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindingAxialRelaxationArm {
    pub relaxing_rigidity_n: f64,
    pub relaxation_time_s: f64,
}

/// Immutable source selection on the same frozen helix and core.
#[derive(Clone, Debug, PartialEq)]
pub struct CloseWoundRelaxation {
    section: CloseWoundSection,
    properties: LayerRelaxationProperties,
    shear_spectrum: GeneralizedMaxwell,
    arms: Vec<WindingAxialRelaxationArm>,
    omega_band_rad_s: (f64, f64),
}

impl CloseWoundSection {
    /// Bind wire-shear history without changing the original equilibrium spring.
    /// Source coefficients must be constant at the held typed environment;
    /// every relaxation time declares its frequency support. An empty arm list
    /// is an explicit elastic selection, requiring an equilibrium source band.
    /// The winding geometry and wavelength restrictions remain unchanged.
    pub fn with_winding_relaxation(
        self,
        properties: LayerRelaxationProperties,
    ) -> Result<CloseWoundRelaxation, MaterialStatePointError> {
        if properties.branches.len() > MAX_SECTION_RELAXATION_BRANCHES {
            return Err(invalid("winding relaxation exceeds the section arm budget"));
        }
        let mut band = self.omega_band_rad_s();
        let mut coefficient = |name: &str, quantity: QuantitySpec, needs_band: bool| {
            let property = self.material.property(name).ok_or_else(|| {
                invalid("winding relaxation needs every selected original shear coefficient")
            })?;
            if property.requirement().quantity() != quantity
                || !property.is_constant_at_fixed_environment()
            {
                return Err(invalid(
                    "winding relaxation requires dimensioned coefficients at the original held environment",
                ));
            }
            match property.angular_frequency_band_rad_s()? {
                Some((lo, hi)) => band = (band.0.max(lo), band.1.min(hi)),
                None if needs_band => {
                    return Err(invalid(
                        "winding relaxation requires explicit source frequency support",
                    ));
                }
                None => (),
            }
            Ok(property.value_si())
        };
        let shear = coefficient(
            self.shear_property(),
            QuantitySpec::dimensional(Pressure::DIMS),
            false,
        )?;
        let equilibrium = coefficient(
            &properties.equilibrium_modulus,
            QuantitySpec::dimensional(Pressure::DIMS),
            properties.branches.is_empty(),
        )?;
        if equilibrium.to_bits() != shear.to_bits() {
            return Err(invalid(
                "winding shear spectrum must retain the original equilibrium modulus",
            ));
        }
        // Reuse the original spring's geometry factor; no second helix model.
        let geometry = positive(
            self.winding_axial_stiffness_n() / shear,
            "positive original helix spring geometry",
        )?;
        let mut terms = Vec::with_capacity(properties.branches.len());
        let mut arms = Vec::with_capacity(properties.branches.len());
        let mut instantaneous = self.axial_stiffness_n();
        for (index, pair) in properties.branches.iter().enumerate() {
            if properties.branches[..index].contains(pair) {
                return Err(invalid("winding relaxation repeats a source-property pair"));
            }
            let delta = coefficient(&pair.0, QuantitySpec::dimensional(Pressure::DIMS), false)?;
            let tau = coefficient(&pair.1, QuantitySpec::dimensional(Time::DIMS), true)?;
            if !delta.is_finite() || delta < 0.0 || !tau.is_finite() || tau <= 0.0 {
                return Err(invalid(
                    "winding relaxation needs nonnegative shear increments and positive source times",
                ));
            }
            let rigidity = delta * geometry;
            instantaneous += rigidity;
            if !rigidity.is_finite()
                || (delta > 0.0 && rigidity == 0.0)
                || !instantaneous.is_finite()
            {
                return Err(invalid("winding relaxation rigidity is unrepresentable"));
            }
            terms.push((delta, tau));
            arms.push(WindingAxialRelaxationArm {
                relaxing_rigidity_n: rigidity,
                relaxation_time_s: tau,
            });
        }
        if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0 {
            return Err(invalid(
                "winding relaxation and original section domains do not intersect",
            ));
        }
        let shear_spectrum = GeneralizedMaxwell::new(equilibrium, terms)
            .map_err(|_| invalid("original winding shear spectrum is not passive and finite"))?;
        Ok(CloseWoundRelaxation {
            section: self,
            properties,
            shear_spectrum,
            arms,
            omega_band_rad_s: band,
        })
    }
}

impl CloseWoundRelaxation {
    pub const fn section(&self) -> &CloseWoundSection {
        &self.section
    }
    pub const fn properties(&self) -> &LayerRelaxationProperties {
        &self.properties
    }
    pub const fn shear_spectrum(&self) -> &GeneralizedMaxwell {
        &self.shear_spectrum
    }
    /// Declaration order, including explicitly zero-strength source arms.
    pub fn arms(&self) -> &[WindingAxialRelaxationArm] {
        &self.arms
    }
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        self.omega_band_rad_s
    }
}
