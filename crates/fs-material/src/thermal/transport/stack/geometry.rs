//! Fixed physical geometry for the original ordered source-backed stack.
//!
//! Radial resistance coordinates and actual inventory metrics are deliberately
//! separate. The finite-volume owner uses equal-volume shells and their volume
//! centroids; changing a conduction coordinate never changes mass or volume.

/// One-dimensional planar or radially symmetric shell geometry. All layers are
/// ordered from the inner face outward. Thicknesses remain physical lengths.
/// Cylinders have insulated endcaps and spheres cover the full solid angle.
/// Positive inner radius is required: a solid centre needs a separate symmetry
/// boundary condition, not a fictitious finite-area inner reservoir.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LayerStackGeometry {
    Planar { area_m2: f64 },
    Cylindrical { inner_radius_m: f64, length_m: f64 },
    Spherical { inner_radius_m: f64 },
}

impl LayerStackGeometry {
    pub const fn is_planar(self) -> bool {
        matches!(self, Self::Planar { .. })
    }

    pub(super) const fn start_m(self) -> f64 {
        match self {
            Self::Planar { .. } => 0.0,
            Self::Cylindrical { inner_radius_m, .. } | Self::Spherical { inner_radius_m } => {
                inner_radius_m
            }
        }
    }

    pub(super) fn reference_area_m2(self) -> Result<f64, &'static str> {
        match self {
            Self::Planar { area_m2 } => checked(area_m2),
            Self::Cylindrical {
                inner_radius_m,
                length_m,
            } => {
                checked(inner_radius_m)?;
                checked(length_m)?;
                self.face_area_m2(inner_radius_m)
            }
            Self::Spherical { inner_radius_m } => {
                checked(inner_radius_m)?;
                self.face_area_m2(inner_radius_m)
            }
        }
    }

    pub(super) fn face_area_m2(self, coordinate_m: f64) -> Result<f64, &'static str> {
        let pi = core::f64::consts::PI;
        checked(match self {
            Self::Planar { area_m2 } => area_m2,
            Self::Cylindrical { length_m, .. } => (2.0 * pi * coordinate_m) * length_m,
            Self::Spherical { .. } => (4.0 * pi * coordinate_m) * coordinate_m,
        })
    }

    /// The two physical face areas, outward-rounded by at most one ulp when
    /// necessary. The native coverage owner subtracts contacts individually;
    /// a nearest-rounded sum must not make a fully covered cell appear negative.
    pub(super) fn cell_surface_area_m2(self, bounds_m: [f64; 2]) -> Result<f64, &'static str> {
        let a = self.face_area_m2(bounds_m[0])?;
        let b = self.face_area_m2(bounds_m[1])?;
        let mut sum = checked(a + b)?;
        if (sum - a) - b < 0.0 || (sum - b) - a < 0.0 {
            sum = checked(sum.next_up())?;
        }
        Ok(sum)
    }

    /// Difference-factored physical shell volume, without nearly equal powers.
    pub(super) fn volume_m3(self, bounds_m: [f64; 2]) -> Result<f64, &'static str> {
        let [a, b] = bounds_m;
        let width = checked(b - a)?;
        let pi = core::f64::consts::PI;
        checked(match self {
            Self::Planar { area_m2 } => area_m2 * width,
            Self::Cylindrical { length_m, .. } => (pi * length_m * width) * (a + b),
            Self::Spherical { .. } => {
                let ratio = a / b;
                ((4.0 * pi / 3.0) * width * b) * b * (1.0 + ratio + ratio * ratio)
            }
        })
    }

    /// A_reference integral(ds / A(s)); this is not a physical layer width.
    pub(super) fn resistance_length_m(self, bounds_m: [f64; 2]) -> Result<f64, &'static str> {
        let [a, b] = bounds_m;
        let width = checked(b - a)?;
        checked(match self {
            Self::Planar { .. } => width,
            Self::Cylindrical { inner_radius_m, .. } => inner_radius_m * (width / a).ln_1p(),
            Self::Spherical { inner_radius_m } => {
                (inner_radius_m / a) * (inner_radius_m / b) * width
            }
        })
    }

    /// Convert an actual interface's areal resistance to the reference area.
    pub(super) fn reference_contact_resistance(
        self,
        resistance: f64,
        interface_m: f64,
    ) -> Result<f64, &'static str> {
        let ratio = self.start_m() / interface_m;
        let factor = match self {
            Self::Planar { .. } => 1.0,
            Self::Cylindrical { .. } => ratio,
            Self::Spherical { .. } => ratio * ratio,
        };
        let value = resistance * factor;
        if value.is_finite() && value >= 0.0 && (resistance == 0.0 || value > 0.0) {
            Ok(value)
        } else {
            Err("radial contact resistance is not representable")
        }
    }

    /// Boundary enclosing the supplied fraction of a layer's physical volume.
    /// Rationalized roots retain thin-shell differences even when a/b is near 1.
    pub(super) fn volume_fraction_position_m(
        self,
        bounds_m: [f64; 2],
        fraction: f64,
    ) -> Result<f64, &'static str> {
        let [a, b] = bounds_m;
        let width = checked(b - a)?;
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err("cell volume fraction must be in [0, 1]");
        }
        if fraction == 0.0 {
            return Ok(a);
        }
        if fraction == 1.0 {
            return Ok(b);
        }
        let ratio = a / b;
        let relative_width = width / b;
        let offset_fraction = match self {
            Self::Planar { .. } => fraction,
            Self::Cylindrical { .. } => {
                let root = (ratio * ratio + fraction * relative_width * (1.0 + ratio)).sqrt();
                fraction * (1.0 + ratio) / (root + ratio)
            }
            Self::Spherical { .. } => {
                let sum = 1.0 + ratio + ratio * ratio;
                let root = (ratio * ratio * ratio + fraction * relative_width * sum).cbrt();
                fraction * sum / (root * root + root * ratio + ratio * ratio)
            }
        };
        let value = a + width * offset_fraction;
        if value.is_finite() && a < value && value < b {
            Ok(value)
        } else {
            Err("radial subdivision is below coordinate resolution")
        }
    }

    /// Volume-weighted radius of a homogeneous reference shell. Difference
    /// forms avoid cancellation and reduce continuously to the planar midpoint.
    pub(super) fn centroid_m(self, bounds_m: [f64; 2]) -> Result<f64, &'static str> {
        let [a, b] = bounds_m;
        let width = checked(b - a)?;
        let ratio = a / b;
        let fraction = match self {
            Self::Planar { .. } => 0.5,
            Self::Cylindrical { .. } => (ratio + 2.0) / (3.0 * (ratio + 1.0)),
            Self::Spherical { .. } => {
                (ratio * ratio + 2.0 * ratio + 3.0) / (4.0 * (ratio * ratio + ratio + 1.0))
            }
        };
        let value = a + width * fraction;
        if value.is_finite() && a < value && value < b {
            Ok(value)
        } else {
            Err("cell centroid is below coordinate resolution")
        }
    }
}

fn checked(value: f64) -> Result<f64, &'static str> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err("layer geometry must be positive, finite and representable")
    }
}
