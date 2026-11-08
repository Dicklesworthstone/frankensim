//! Shared reference-frame kinetic integral for prismatic 3-D beam segments.

use super::{BeamSection, bad};
use crate::PlateError;

mod section;
pub use section::BeamSectionInertia;

#[derive(Clone, Copy)]
pub(in crate::shell) enum BeamReferenceInertia {
    Density([f64; 2]),
    Section([BeamSectionInertia; 2]),
}

impl BeamReferenceInertia {
    pub(in crate::shell) fn validate(self) -> Result<(), PlateError> {
        match self {
            Self::Density(rho) => validate_reference_density(rho),
            Self::Section(_) => Ok(()), // Immutable, admitted physical moments.
        }
    }

    pub(in crate::shell) fn local_mass(
        self,
        section: BeamSection,
        length: f64,
        phi: [f64; 2],
        mass: BeamMass,
    ) -> Result<[f64; 144], PlateError> {
        match self {
            Self::Density(rho) => {
                local_mass_with_reference_density(section, length, phi, mass, rho)
            }
            Self::Section(moments) => section::local_mass(section, length, phi, mass, moments),
        }
    }
}

#[derive(Clone, Copy)]
pub(in crate::shell) enum BeamReferenceInertiaField<'a> {
    Density(&'a [[f64; 2]]),
    Section(&'a [[BeamSectionInertia; 2]]),
}

impl BeamReferenceInertiaField<'_> {
    pub(in crate::shell) fn validate(self, segments: usize) -> Result<(), PlateError> {
        let len = match self {
            Self::Density(v) => v.len(),
            Self::Section(v) => v.len(),
        };
        if len != segments {
            return Err(bad(
                "beam reference inertia must cover every original segment",
            ));
        }
        if let Self::Density(values) = self {
            for &value in values {
                validate_reference_density(value)?;
            }
        }
        Ok(())
    }

    pub(in crate::shell) fn at(self, index: usize) -> BeamReferenceInertia {
        match self {
            Self::Density(v) => BeamReferenceInertia::Density(v[index]),
            Self::Section(v) => BeamReferenceInertia::Section(v[index]),
        }
    }
}

// Shared original four-point Gauss rule; P1 section moments increase the
// degree-six velocity product to degree seven without changing the sites.
const GAUSS: [(f64, f64); 4] = [
    (0.06943184420297371, 0.17392742256872693),
    (0.33000947820757187, 0.32607257743127307),
    (0.6699905217924281, 0.32607257743127307),
    (0.9305681557970262, 0.17392742256872693),
];

// Original independent static Timoshenko velocity and section-rotation fields.
fn bending_fields(l: f64, p: f64, t: f64) -> ([f64; 4], [f64; 4]) {
    let r = 1.0 / (1.0 + p);
    let f = r * (3.0 * t * t - 2.0 * t * t * t) + (p * r) * t;
    let a = 6.0 * r * (t - t * t) / l;
    let n = [
        1.0 - f,
        l * (t - 0.5 * t * t - 0.5 * f),
        f,
        l * (0.5 * t * t - 0.5 * f),
    ];
    let h = [-a, 1.0 - t - 0.5 * l * a, a, t - 0.5 * l * a];
    (n, h)
}

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
    local_mass_impl(s, l, phi, mass, None)
}

pub(super) fn validate_reference_density(density: [f64; 2]) -> Result<(), PlateError> {
    if density.iter().any(|rho| !rho.is_finite() || *rho <= 0.0) {
        return Err(bad(
            "beam reference density must be finite and positive at both endpoints",
        ));
    }
    Ok(())
}

/// P1 original density on the same prismatic geometry and kinetic fields.
/// A uniform field delegates to the scalar owner, preserving its arithmetic.
pub(super) fn local_mass_with_reference_density(
    s: BeamSection,
    l: f64,
    phi: [f64; 2],
    mass: BeamMass,
    reference_density_kg_m3: [f64; 2],
) -> Result<[f64; 144], PlateError> {
    validate_reference_density(reference_density_kg_m3)?;
    if reference_density_kg_m3[0] == reference_density_kg_m3[1] {
        return local_mass(
            BeamSection {
                density_kg_m3: reference_density_kg_m3[0],
                ..s
            },
            l,
            phi,
            mass,
        );
    }
    local_mass_impl(s, l, phi, mass, Some(reference_density_kg_m3))
}

fn local_mass_impl(
    s: BeamSection,
    l: f64,
    phi: [f64; 2],
    mass: BeamMass,
    reference_density: Option<[f64; 2]>,
) -> Result<[f64; 144], PlateError> {
    let mut m = [0.0; 144];
    if mass == BeamMass::Lumped {
        // Preserve the arithmetic of the original bonded ShellBeam pencil.
        let half = s.density_kg_m3 * l / 2.0;
        // Integral rho*N_i dx preserves both total inventory and its first
        // axial moment. Equal endpoint halves would lose the density centroid.
        let partition = reference_density.map_or([half; 2], |rho| {
            [
                l * (rho[0] / 3.0 + rho[1] / 6.0),
                l * (rho[0] / 6.0 + rho[1] / 3.0),
            ]
        });
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
                m[12 * i + i] = partition[node] * inertia;
            }
        }
    } else {
        if phi.iter().any(|p| !p.is_finite() || *p < 0.0) {
            return Err(bad("unresolved beam shear-to-bending ratio"));
        }
        // Axial translation and torsion have linear fields. Torsional inertia is
        // the polar area moment Iy+Iz, never the Saint-Venant stiffness constant J.
        for (a, b, inertia) in [(0, 6, s.area_m2), (3, 9, s.iy_m4 + s.iz_m4)] {
            let (left, off, right) = if let Some(rho) = reference_density {
                // Exact P1-density/P1-velocity Gram entries, without replacing
                // the density by its mean or forming an overflowing sum first.
                (
                    l * inertia * (rho[0] / 4.0 + rho[1] / 12.0),
                    l * inertia * (rho[0] / 12.0 + rho[1] / 12.0),
                    l * inertia * (rho[0] / 12.0 + rho[1] / 4.0),
                )
            } else {
                let off = s.density_kg_m3 * l * inertia / 6.0;
                (2.0 * off, off, 2.0 * off)
            };
            m[12 * a + a] = left;
            m[12 * b + b] = right;
            m[12 * a + b] = off;
            m[12 * b + a] = off;
        }
        // Four Gauss points integrate the degree-six velocity products, or
        // degree seven with P1 density, without any additional quadrature.
        for (slots, signs, inertia, p) in [
            ([1, 5, 7, 11], [1.0; 4], s.iz_m4, phi[0]),
            ([2, 4, 8, 10], [1.0, -1.0, 1.0, -1.0], s.iy_m4, phi[1]),
        ] {
            for (t, weight) in GAUSS {
                let rho = reference_density
                    .map_or(s.density_kg_m3, |rho| (1.0 - t) * rho[0] + t * rho[1]);
                // Exact static Timoshenko fields: V=GA_s(v'-theta) is constant,
                // EI theta''=-V. Cubic v and quadratic theta interpolate the four
                // endpoint DOFs; at phi=0, theta=v' and v is cubic Hermite.
                let (n, h) = bending_fields(l, p, t);
                for i in 0..4 {
                    for j in i..4 {
                        let value = rho
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
