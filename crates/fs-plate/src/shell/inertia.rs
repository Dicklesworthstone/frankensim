//! Kinetic interpolation for the existing flat-facet CST/DKT shell owner.
//!
//! P1 midsurface velocity and P1 physical director rotation are integrated
//! through the section thickness. This is not a new DKT displacement field.

use super::FacetGeometry;
use crate::{PlateError, PlateSection};

/// Mass discretization of the flat-facet CST/DKT shell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CstDktMass {
    /// Original diagonal translation/rotation mass, retained for frozen pencils.
    Lumped,
    /// Exact P1 midsurface/director kinetic integral and a relative drilling gauge.
    /// The gauge retains the original `rho*h^3/12` scale but vanishes for rigid
    /// motion; it is numerical stabilization, not material normal-spin inertia.
    Consistent,
}

/// Local `[u,v,w,theta_x,theta_y,theta_z]` coordinates at each vertex.
pub(super) fn consistent_mass(
    g: &FacetGeometry,
    s: &PlateSection,
) -> Result<[f64; 324], PlateError> {
    let total = s.density * s.thickness * g.area_m2;
    let rotary = total * s.thickness * s.thickness / 12.0;
    if !total.is_finite() || total <= 0.0 || !rotary.is_finite() || rotary <= 0.0 {
        return Err(PlateError::BadSection {
            what: "consistent shell inertia is unrepresentable",
        });
    }
    let mut m = [0.0; 324];
    // Integral Ni*Nj / area = (1 + delta_ij)/12. Rotation about the
    // facet normal does not move the through-thickness material director.
    for i in 0..3 {
        for j in 0..3 {
            let gram = if i == j { 1.0 / 6.0 } else { 1.0 / 12.0 };
            for c in 0..5 {
                m[18 * (6 * i + c) + 6 * j + c] = gram * if c < 3 { total } else { rotary };
            }
        }
    }
    // The original drilling stiffness penalizes theta_z - (v_x-u_y)/2.
    // Stabilize its RATE, so an actual rigid normal spin acquires no invented
    // section inertia. Positive degree-two edge quadrature is exact for P1.
    let slots = [0, 1, 5, 6, 7, 11, 12, 13, 17];
    for h in [[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]] {
        let spin: [f64; 9] = std::array::from_fn(|i| match i % 3 {
            0 => 0.5 * g.gradient[i / 3][1],
            1 => -0.5 * g.gradient[i / 3][0],
            _ => h[i / 3],
        });
        for i in 0..9 {
            for j in i..9 {
                let value = rotary / 3.0 * spin[i] * spin[j];
                m[18 * slots[i] + slots[j]] += value;
                if i != j {
                    m[18 * slots[j] + slots[i]] += value;
                }
            }
        }
    }
    if m.iter().any(|v| !v.is_finite()) {
        return Err(PlateError::BadSection {
            what: "consistent shell kinetic integral overflow",
        });
    }
    Ok(m)
}
