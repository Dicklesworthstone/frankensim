//! Shared centroidal kinetic integral for prismatic 3-D beam segments.

use super::{BeamSection, bad};
use crate::PlateError;

/// Inertia discretization of a straight beam; both choices include section rotary inertia.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BeamMass {
    /// Integrate the beam's displacement and independent section-rotation fields exactly.
    /// This is the normal segment/chain path, including beam-to-shell ties.
    Consistent,
    /// Endpoint quadrature retained for the frozen bonded piano-board pencil.
    /// It preserves total mass but overestimates a single segment's rigid transverse
    /// orbital inertia by a factor of three. Refine the mesh or use `Consistent`.
    Lumped,
}

/// `phi` is [12 E Iz/(G As_y L²), 12 E Iy/(G As_z L²)]; zero gives Rayleigh/EB inertia.
pub(super) fn local_mass(
    s: BeamSection,
    l: f64,
    phi: [f64; 2],
    mass: BeamMass,
) -> Result<[f64; 144], PlateError> {
    let mut m = [0.0; 144];
    if mass == BeamMass::Lumped {
        // Preserve the arithmetic of the original bonded ShellBeam pencil.
        let half = s.density_kg_m3 * l / 2.0;
        for node in 0..2 {
            for (dof, inertia) in [
                s.area_m2,
                s.area_m2,
                s.area_m2,
                s.iy_m4 + s.iz_m4,
                s.iy_m4,
                s.iz_m4,
            ]
            .into_iter()
            .enumerate()
            {
                let i = 6 * node + dof;
                m[12 * i + i] = half * inertia;
            }
        }
    } else {
        if phi.iter().any(|p| !p.is_finite() || *p < 0.0) {
            return Err(bad("unresolved beam shear-to-bending ratio"));
        }
        // Axial translation and torsion have linear fields. Torsional inertia is
        // the polar area moment Iy+Iz, never the Saint-Venant stiffness constant J.
        for (a, b, inertia) in [(0, 6, s.area_m2), (3, 9, s.iy_m4 + s.iz_m4)] {
            let off = s.density_kg_m3 * l * inertia / 6.0;
            m[12 * a + a] = 2.0 * off;
            m[12 * b + b] = 2.0 * off;
            m[12 * a + b] = off;
            m[12 * b + a] = off;
        }
        // Four Gauss points integrate the degree-six kinetic polynomial exactly.
        const GAUSS: [(f64, f64); 4] = [
            (0.06943184420297371, 0.17392742256872693),
            (0.33000947820757187, 0.32607257743127307),
            (0.6699905217924281, 0.32607257743127307),
            (0.9305681557970262, 0.17392742256872693),
        ];
        for (slots, signs, inertia, p) in [
            ([1, 5, 7, 11], [1.0; 4], s.iz_m4, phi[0]),
            ([2, 4, 8, 10], [1.0, -1.0, 1.0, -1.0], s.iy_m4, phi[1]),
        ] {
            let r = 1.0 / (1.0 + p);
            for (t, weight) in GAUSS {
                // Exact static Timoshenko fields: V=GA_s(v'-theta) is constant,
                // EI theta''=-V. Cubic v and quadratic theta interpolate the four
                // endpoint DOFs; at phi=0, theta=v' and v is cubic Hermite.
                let f = r * (3.0 * t * t - 2.0 * t * t * t) + (p * r) * t;
                let a = 6.0 * r * (t - t * t) / l;
                let n = [
                    1.0 - f,
                    l * (t - 0.5 * t * t - 0.5 * f),
                    f,
                    l * (0.5 * t * t - 0.5 * f),
                ];
                let h = [-a, 1.0 - t - 0.5 * l * a, a, t - 0.5 * l * a];
                for i in 0..4 {
                    for j in i..4 {
                        let value = s.density_kg_m3
                            * l
                            * weight
                            * (s.area_m2 * n[i] * n[j] + inertia * h[i] * h[j])
                            * signs[i]
                            * signs[j];
                        let (row, col) = (slots[i], slots[j]);
                        m[12 * row + col] += value;
                        if row != col {
                            m[12 * col + row] += value;
                        }
                    }
                }
            }
        }
    }
    if m.iter().any(|x| !x.is_finite()) {
        return Err(bad("beam inertia overflow"));
    }
    Ok(m)
}
