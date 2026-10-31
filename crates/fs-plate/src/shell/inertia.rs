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

/// The same kinetic field with thickness interpolated from the measured vertices.
/// The section supplies density; its facet-average thickness still owns stiffness.
/// Positive degree-five triangle quadrature integrates h*Ni*Nj and h^3*Ni*Nj
/// exactly (up to roundoff), including the relative drilling rate. The material
/// occupies -h(x,y)/2 <= z <= h(x,y)/2 about this flat reference facet.
pub(super) fn consistent_mass_with_thickness(
    g: &FacetGeometry,
    s: &PlateSection,
    thickness: [f64; 3],
) -> Result<[f64; 324], PlateError> {
    if thickness[0] == thickness[1] && thickness[1] == thickness[2] {
        // Preserve the original section-average arithmetic, including rounding,
        // for every homogeneous facet and its existing modal pencils.
        return consistent_mass(g, s);
    }
    // Symmetric seven-point degree-five rule, weights normalized to unit area.
    // The two orbits have repeated coordinate (6 +/- sqrt(15))/21 and
    // weight (155 +/- sqrt(15))/1200. All points and weights are positive.
    const A: f64 = 0.4701420641051151;
    const B: f64 = 0.05971587178976982;
    const C: f64 = 0.10128650732345634;
    const D: f64 = 0.7974269853530873;
    const W: f64 = 0.1323941527885062;
    const V: f64 = 0.12593918054482714;
    let rule = [
        ([1.0 / 3.0; 3], 0.225),
        ([B, A, A], W),
        ([A, B, A], W),
        ([A, A, B], W),
        ([D, C, C], V),
        ([C, D, C], V),
        ([C, C, D], V),
    ];
    let mut m = [0.0; 324];
    let slots = [0, 1, 5, 6, 7, 11, 12, 13, 17];
    for (n, weight) in rule {
        let h = n[0] * thickness[0] + n[1] * thickness[1] + n[2] * thickness[2];
        let translation = s.density * h * g.area_m2 * weight;
        let rotary = translation * h * h / 12.0;
        if !translation.is_finite() || translation <= 0.0 || !rotary.is_finite() || rotary <= 0.0 {
            return Err(PlateError::BadSection {
                what: "variable-thickness shell inertia is unrepresentable",
            });
        }
        for i in 0..3 {
            for j in i..3 {
                for c in 0..5 {
                    let value = n[i] * n[j] * if c < 3 { translation } else { rotary };
                    m[18 * (6 * i + c) + 6 * j + c] += value;
                    if i != j {
                        m[18 * (6 * j + c) + 6 * i + c] += value;
                    }
                }
            }
        }
        let spin: [f64; 9] = std::array::from_fn(|i| match i % 3 {
            0 => 0.5 * g.gradient[i / 3][1],
            1 => -0.5 * g.gradient[i / 3][0],
            _ => n[i / 3],
        });
        for i in 0..9 {
            for j in i..9 {
                let value = rotary * spin[i] * spin[j];
                m[18 * slots[i] + slots[j]] += value;
                if i != j {
                    m[18 * slots[j] + slots[i]] += value;
                }
            }
        }
    }
    if m.iter().any(|v| !v.is_finite()) {
        return Err(PlateError::BadSection {
            what: "variable-thickness shell kinetic integral overflow",
        });
    }
    Ok(m)
}
