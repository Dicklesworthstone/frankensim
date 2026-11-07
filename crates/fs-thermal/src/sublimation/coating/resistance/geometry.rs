//! Conductive coating geometry derived from the original conserved phase mass.
//!
//! A fixed substrate supports a uniform planar, spherical or cylindrical layer.
//! [`RetreatGeometry`] remains the only mass-to-thickness/area owner. The layer's
//! full-thickness, constant-conductivity resistance and an authored contact
//! resistance supply the original two-temperature surface solver with one
//! conductance. The support contact area is distinct from the exposed outer
//! area on curved layers.
//!
//! These are the steady Fourier resistances `integral dr/(k*A(r))`; see
//! [MIT Unified Thermodynamics, section 16.5](https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node119.html).
//! Using that resistance with one stored deposit temperature is a two-node
//! transient approximation, not a resolved through-layer temperature field.
//! A caller refreshes the geometry and conductivity between time intervals.
//! No thermal expansion, porosity, densification, nucleation or gas-volume work
//! is inferred. Zero coating mass requires the separate interface-transition
//! owner; zero conductivity is not used to invent a material law.

use crate::phase_geometry::{FrontGeometry, RetreatGeometry};
use crate::sublimation::SublimationError;
use core::f64::consts::PI;

#[derive(Clone, Copy, Debug, PartialEq)]
enum Shape {
    Planar,
    Spherical { core_radius_m: f64 },
    Cylindrical { core_radius_m: f64, length_m: f64 },
}

/// Constant-density supported geometry, with fixed inner contact area.
///
/// The original [`RetreatGeometry`] derives every moving dimension from deposit
/// mass. Cylinder end caps exchange neither mass nor heat in this reduction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductiveCoatingGeometry {
    geometry: RetreatGeometry,
    shape: Shape,
    contact_area_m2: f64,
}

/// One mass-owned geometry and its full-thickness thermal circuit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductiveCoatingEndpoint {
    /// Original coating-only mass/volume, outer area and thickness.
    pub front: FrontGeometry,
    /// Fixed support/coating contact area; not the growing outer area.
    pub contact_area_m2: f64,
    /// Constant-conductivity layer resistance [K/W], excluding contact.
    pub layer_resistance_k_w: f64,
    /// Authored areal contact resistance divided by the inner area [K/W].
    pub contact_resistance_k_w: f64,
    /// `1 / (layer_resistance_k_w + contact_resistance_k_w)` [W/K].
    pub conductance_w_k: f64,
}

fn invalid(message: &'static str) -> SublimationError {
    SublimationError::Invalid(message)
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

impl ConductiveCoatingGeometry {
    /// Uniform layer on one fixed footprint, with inert sides and back face.
    pub fn planar(density_kg_m3: f64, footprint_m2: f64) -> Result<Self, SublimationError> {
        Ok(Self {
            geometry: RetreatGeometry::layer(density_kg_m3, footprint_m2)?,
            shape: Shape::Planar,
            contact_area_m2: footprint_m2,
        })
    }

    /// Uniform shell around a fixed spherical substrate.
    pub fn spherical(density_kg_m3: f64, core_radius_m: f64) -> Result<Self, SublimationError> {
        let geometry = RetreatGeometry::spherical_coating(density_kg_m3, core_radius_m)?;
        let contact_area_m2 = (4.0 * PI * core_radius_m) * core_radius_m;
        if !positive(contact_area_m2) {
            return Err(invalid("spherical coating contact area is unrepresentable"));
        }
        Ok(Self {
            geometry,
            shape: Shape::Spherical { core_radius_m },
            contact_area_m2,
        })
    }

    /// Radial coating of a fixed cylinder; end caps are thermally and materially inert.
    pub fn cylindrical(
        density_kg_m3: f64,
        core_radius_m: f64,
        length_m: f64,
    ) -> Result<Self, SublimationError> {
        let geometry =
            RetreatGeometry::cylindrical_coating(density_kg_m3, core_radius_m, length_m)?;
        let contact_area_m2 = (2.0 * PI * core_radius_m) * length_m;
        if !positive(contact_area_m2) {
            return Err(invalid(
                "cylindrical coating contact area is unrepresentable",
            ));
        }
        Ok(Self {
            geometry,
            shape: Shape::Cylindrical {
                core_radius_m,
                length_m,
            },
            contact_area_m2,
        })
    }

    /// Original mass/volume/recession owner, including substrate displacement.
    #[must_use]
    pub const fn geometry(self) -> RetreatGeometry {
        self.geometry
    }

    #[must_use]
    pub const fn contact_area_m2(self) -> f64 {
        self.contact_area_m2
    }

    /// Keep the phase law's density consistent with this constant-density geometry.
    pub fn admit_density(self, density_kg_m3: f64) -> Result<(), SublimationError> {
        self.geometry.admit_density(density_kg_m3)
    }

    /// Refresh outer area, thickness and conductance from the same actual mass.
    ///
    /// `conductivity_w_m_k` is a positive original material property frozen over
    /// this interval. `contact_resistance_m2_k_w` is a finite nonnegative property
    /// of the support/coating interface. It is referred to the INNER contact
    /// area and is added in series with the layer resistance. No gas-film
    /// resistance is included. Invalid or unrepresentable circuits refuse;
    /// they are not clipped into perfect contact or thermal insulation.
    ///
    /// Spherical resistance uses the factored thickness rather than subtracting
    /// reciprocal radii. The cylindrical logarithm compensates the rounding of
    /// `1 + thickness/core_radius`, retaining films thinner than one radius ULP.
    pub fn at_mass(
        self,
        mass_kg: f64,
        conductivity_w_m_k: f64,
        contact_resistance_m2_k_w: f64,
    ) -> Result<ConductiveCoatingEndpoint, SublimationError> {
        if !positive(mass_kg)
            || !positive(conductivity_w_m_k)
            || !contact_resistance_m2_k_w.is_finite()
            || contact_resistance_m2_k_w < 0.0
        {
            return Err(invalid(
                "positive coating mass and conductivity, nonnegative contact resistance required",
            ));
        }
        let front = self.geometry.at_mass(mass_kg)?;
        let thickness = front.extent_m;
        let layer_resistance_k_w = match self.shape {
            Shape::Planar => (thickness / self.contact_area_m2) / conductivity_w_m_k,
            Shape::Spherical { core_radius_m } => {
                let outer_radius = core_radius_m + thickness;
                // (1/ri - 1/ro)/(4*pi*k), rationalized before evaluation.
                ((thickness / outer_radius) / core_radius_m / (4.0 * PI)) / conductivity_w_m_k
            }
            Shape::Cylindrical {
                core_radius_m,
                length_m,
            } => {
                let ratio = thickness / core_radius_m;
                let rounded = 1.0 + ratio;
                let logarithm = if rounded.is_finite() {
                    // Original deterministic logarithm plus the first-order
                    // correction for the rounding of 1+x. For 1+x==1 this
                    // reduces to x, without losing the physical layer.
                    fs_math::det::ln(rounded) + (ratio - (rounded - 1.0)) / rounded
                } else {
                    // Widely separated radii cannot suffer cancellation here.
                    fs_math::det::ln(core_radius_m + thickness) - fs_math::det::ln(core_radius_m)
                };
                (logarithm / length_m / (2.0 * PI)) / conductivity_w_m_k
            }
        };
        let contact_resistance_k_w = contact_resistance_m2_k_w / self.contact_area_m2;
        let total_resistance = layer_resistance_k_w + contact_resistance_k_w;
        let conductance_w_k = 1.0 / total_resistance;
        if !positive(layer_resistance_k_w)
            || !contact_resistance_k_w.is_finite()
            || contact_resistance_m2_k_w > 0.0 && contact_resistance_k_w == 0.0
            || !positive(total_resistance)
            || !positive(conductance_w_k)
        {
            return Err(invalid("coating thermal circuit is unrepresentable"));
        }
        Ok(ConductiveCoatingEndpoint {
            front,
            contact_area_m2: self.contact_area_m2,
            layer_resistance_k_w,
            contact_resistance_k_w,
            conductance_w_k,
        })
    }
}

#[cfg(test)]
#[path = "geometry_tests.rs"]
mod tests;
