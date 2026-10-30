//! Exact resistance coordinates for radial, source-free Fourier conduction.
//!
//! With area A(r), steady power satisfies integral(k dT) = Q integral(dr/A).
//! Multiplying each geometric integral by the actual contact area produces
//! equivalent planar distances. This changes only geometry: the original
//! conductivity integrals, contact solve and endpoint derivatives remain owners.

use super::{ConductivityCurve, DepositionError, TwoLayerConduction};

/// Admitted planar, cylindrical or spherical conduction between two temperature
/// traces and their intervening material interface. The trace temperatures may
/// be finite-cell unknowns or prescribed values; this geometry stores no heat.
///
/// Radial traces must straddle the positive interface radius. Either a/b order
/// is allowed. Cylinders have insulated ends, spheres have full solid angle;
/// scalar k must describe the radial direction. Radial symmetry, fixed geometry
/// and source-free quasi-steady legs are explicit model assumptions. This is
/// not a convection law, fluid motion, a point source at r=0, or an assertion
/// that unresolved layer storage is negligible in every transient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerGeometry {
    equivalent_distances_m: [f64; 2],
    contact_area_m2: f64,
    coordinates: Coordinates,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Coordinates {
    Planar,
    Cylindrical {
        temperature_radii_m: [f64; 2],
        contact_radius_m: f64,
        length_m: f64,
    },
    Spherical {
        temperature_radii_m: [f64; 2],
        contact_radius_m: f64,
    },
}

impl LayerGeometry {
    /// Homogeneous planar legs of the stated physical lengths and common area.
    pub fn planar(distances_m: [f64; 2], area_m2: f64) -> Result<Self, DepositionError> {
        Self::admit(distances_m, area_m2, Coordinates::Planar)
    }

    /// Exact annular-cylinder metrics: d_i = r_contact |ln(r_i/r_contact)|.
    /// Total power uses the actual cylindrical contact area 2*pi*r_contact*L.
    /// `ln_1p` retains close-radius differences in thin coatings.
    pub fn cylindrical(
        temperature_radii_m: [f64; 2],
        contact_radius_m: f64,
        length_m: f64,
    ) -> Result<Self, DepositionError> {
        admit_radii(temperature_radii_m, contact_radius_m)?;
        if !positive(length_m) {
            return Err(invalid(
                "cylindrical conduction needs positive finite length",
            ));
        }
        let distances = temperature_radii_m.map(|r| {
            let low = r.min(contact_radius_m);
            let high = r.max(contact_radius_m);
            contact_radius_m * ((high - low) / low).ln_1p()
        });
        let area = (2.0 * core::f64::consts::PI * contact_radius_m) * length_m;
        Self::admit(
            distances,
            area,
            Coordinates::Cylindrical {
                temperature_radii_m,
                contact_radius_m,
                length_m,
            },
        )
    }

    /// Exact spherical metrics: d_i = r_contact*|r_contact-r_i|/r_i.
    /// Total power uses the actual spherical contact area 4*pi*r_contact².
    /// The difference form avoids subtracting nearly equal inverse radii.
    pub fn spherical(
        temperature_radii_m: [f64; 2],
        contact_radius_m: f64,
    ) -> Result<Self, DepositionError> {
        admit_radii(temperature_radii_m, contact_radius_m)?;
        let distances =
            temperature_radii_m.map(|r| contact_radius_m * ((contact_radius_m - r).abs() / r));
        let area = (4.0 * core::f64::consts::PI * contact_radius_m) * contact_radius_m;
        Self::admit(
            distances,
            area,
            Coordinates::Spherical {
                temperature_radii_m,
                contact_radius_m,
            },
        )
    }

    fn admit(
        equivalent_distances_m: [f64; 2],
        contact_area_m2: f64,
        coordinates: Coordinates,
    ) -> Result<Self, DepositionError> {
        if !equivalent_distances_m.into_iter().all(positive) || !positive(contact_area_m2) {
            return Err(invalid(
                "conduction geometry is nonpositive or unrepresentable",
            ));
        }
        Ok(Self {
            equivalent_distances_m,
            contact_area_m2,
            coordinates,
        })
    }

    /// Physical material-interface area, also used for its areal resistance.
    pub const fn contact_area_m2(self) -> f64 {
        self.contact_area_m2
    }

    /// A_contact integral(ds/A(s)) on each ordered leg [m]. These are physical
    /// lengths only for a planar path. Do not use them for mass or cell volume.
    pub const fn equivalent_distances_m(self) -> [f64; 2] {
        self.equivalent_distances_m
    }

    /// Original a/b radial temperature coordinates and interface radius, absent
    /// for a planar path. The caller retains the finite-cell volume discretization.
    pub const fn radial_coordinates_m(self) -> Option<([f64; 2], f64)> {
        match self.coordinates {
            Coordinates::Planar => None,
            Coordinates::Cylindrical {
                temperature_radii_m,
                contact_radius_m,
                ..
            }
            | Coordinates::Spherical {
                temperature_radii_m,
                contact_radius_m,
            } => Some((temperature_radii_m, contact_radius_m)),
        }
    }

    /// Bind the original k(T) laws and area-specific interface resistance to the
    /// shared nonlinear flux owner. Successful construction allocates nothing.
    pub fn conduction<'a>(
        self,
        curves: [&'a ConductivityCurve; 2],
        contact_resistance_m2_k_w: f64,
    ) -> Result<TwoLayerConduction<'a>, DepositionError> {
        TwoLayerConduction::new(
            curves,
            self.equivalent_distances_m,
            self.contact_area_m2,
            contact_resistance_m2_k_w,
        )
    }
}

fn admit_radii(radii: [f64; 2], contact: f64) -> Result<(), DepositionError> {
    if !radii.into_iter().all(positive)
        || !positive(contact)
        || !((radii[0] < contact && contact < radii[1])
            || (radii[1] < contact && contact < radii[0]))
    {
        return Err(invalid(
            "radial temperature traces must straddle a positive contact radius",
        ));
    }
    Ok(())
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn invalid(reason: &'static str) -> DepositionError {
    DepositionError::InvalidInput(reason)
}
